% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRP1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:20:22
	% EndTime: 2022-01-20 10:20:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:20:23
	% EndTime: 2022-01-20 10:20:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:20:22
	% EndTime: 2022-01-20 10:20:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t43 = pkin(1) * qJD(1);
	t40 = qJ(1) + qJ(2);
	t37 = sin(t40);
	t38 = cos(t40);
	t39 = qJD(1) + qJD(2);
	t42 = (-r_i_i_C(1) * t38 + r_i_i_C(2) * t37) * t39;
	t41 = (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * t39;
	t1 = [-cos(qJ(1)) * t43 + t42, t42, 0, 0, 0; -sin(qJ(1)) * t43 + t41, t41, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:20:22
	% EndTime: 2022-01-20 10:20:23
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (42->9), mult. (28->11), div. (0->0), fcn. (14->6), ass. (0->9)
	t49 = pkin(1) * qJD(1);
	t46 = qJ(1) + qJ(2);
	t42 = pkin(8) + t46;
	t40 = sin(t42);
	t41 = cos(t42);
	t45 = qJD(1) + qJD(2);
	t48 = (-pkin(2) * cos(t46) - r_i_i_C(1) * t41 + t40 * r_i_i_C(2)) * t45;
	t47 = (-pkin(2) * sin(t46) - r_i_i_C(1) * t40 - r_i_i_C(2) * t41) * t45;
	t1 = [-cos(qJ(1)) * t49 + t48, t48, 0, 0, 0; -sin(qJ(1)) * t49 + t47, t47, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:20:23
	% EndTime: 2022-01-20 10:20:23
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (141->22), mult. (112->35), div. (0->0), fcn. (68->8), ass. (0->21)
	t66 = r_i_i_C(3) + pkin(7);
	t48 = qJ(1) + qJ(2);
	t44 = pkin(8) + t48;
	t42 = sin(t44);
	t43 = cos(t44);
	t50 = cos(qJ(4));
	t59 = qJD(4) * t50;
	t47 = qJD(1) + qJD(2);
	t49 = sin(qJ(4));
	t63 = t47 * t49;
	t65 = t42 * t59 + t43 * t63;
	t62 = t47 * t50;
	t61 = pkin(1) * qJD(1);
	t60 = qJD(4) * t49;
	t58 = t42 * t63;
	t56 = t42 * t60;
	t54 = -r_i_i_C(1) * t50 - pkin(3);
	t53 = (-r_i_i_C(1) * t49 - r_i_i_C(2) * t50) * qJD(4);
	t52 = r_i_i_C(1) * t56 + (-pkin(2) * cos(t48) + t54 * t43 - t66 * t42) * t47 + t65 * r_i_i_C(2);
	t51 = r_i_i_C(2) * t58 + (-pkin(2) * sin(t48) + t54 * t42) * t47 + (t66 * t47 + t53) * t43;
	t1 = [-cos(qJ(1)) * t61 + t52, t52, 0, (t42 * t62 + t43 * t60) * r_i_i_C(2) + (-t43 * t59 + t58) * r_i_i_C(1), 0; -sin(qJ(1)) * t61 + t51, t51, 0, (-t43 * t62 + t56) * r_i_i_C(2) - t65 * r_i_i_C(1), 0; 0, 0, 0, t53, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:20:23
	% EndTime: 2022-01-20 10:20:23
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (194->21), mult. (148->26), div. (0->0), fcn. (95->8), ass. (0->21)
	t54 = qJD(1) + qJD(2);
	t80 = pkin(2) * t54;
	t57 = sin(qJ(4));
	t58 = cos(qJ(4));
	t74 = pkin(4) + r_i_i_C(1);
	t75 = r_i_i_C(2) * t58 + t74 * t57;
	t61 = t75 * qJD(4);
	t79 = -t61 - t54 * (-qJ(5) - pkin(7));
	t78 = r_i_i_C(2) * t57 - t74 * t58;
	t76 = qJD(5) + (-pkin(3) + t78) * t54;
	t55 = qJ(1) + qJ(2);
	t51 = pkin(8) + t55;
	t48 = sin(t51);
	t71 = t54 * t48;
	t49 = cos(t51);
	t70 = t54 * t49;
	t68 = pkin(1) * qJD(1);
	t62 = qJD(4) * t78;
	t60 = -cos(t55) * t80 + t76 * t49 + (-r_i_i_C(3) * t54 - t79) * t48;
	t59 = r_i_i_C(3) * t70 - sin(t55) * t80 + t79 * t49 + t76 * t48;
	t1 = [-cos(qJ(1)) * t68 + t60, t60, 0, t49 * t62 + t75 * t71, t70; -sin(qJ(1)) * t68 + t59, t59, 0, t48 * t62 - t70 * t75, t71; 0, 0, 0, -t61, 0;];
	JaD_transl = t1;
end