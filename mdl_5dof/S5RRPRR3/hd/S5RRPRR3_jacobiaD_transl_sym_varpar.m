% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:34:48
	% EndTime: 2022-01-20 10:34:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:34:48
	% EndTime: 2022-01-20 10:34:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:34:48
	% EndTime: 2022-01-20 10:34:48
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
	% StartTime: 2022-01-20 10:34:48
	% EndTime: 2022-01-20 10:34:48
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (42->9), mult. (28->11), div. (0->0), fcn. (14->6), ass. (0->9)
	t49 = pkin(1) * qJD(1);
	t46 = qJ(1) + qJ(2);
	t42 = pkin(9) + t46;
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
	% StartTime: 2022-01-20 10:34:48
	% EndTime: 2022-01-20 10:34:48
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (96->13), mult. (44->15), div. (0->0), fcn. (22->8), ass. (0->13)
	t60 = pkin(1) * qJD(1);
	t55 = qJ(1) + qJ(2);
	t54 = qJD(1) + qJD(2);
	t51 = pkin(9) + t55;
	t49 = qJ(4) + t51;
	t45 = sin(t49);
	t46 = cos(t49);
	t50 = qJD(4) + t54;
	t59 = (-r_i_i_C(1) * t46 + r_i_i_C(2) * t45) * t50;
	t58 = (-r_i_i_C(1) * t45 - r_i_i_C(2) * t46) * t50;
	t57 = (-pkin(2) * cos(t55) - pkin(3) * cos(t51)) * t54 + t59;
	t56 = (-pkin(2) * sin(t55) - pkin(3) * sin(t51)) * t54 + t58;
	t1 = [-cos(qJ(1)) * t60 + t57, t57, 0, t59, 0; -sin(qJ(1)) * t60 + t56, t56, 0, t58, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:34:48
	% EndTime: 2022-01-20 10:34:48
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (291->26), mult. (156->39), div. (0->0), fcn. (94->10), ass. (0->24)
	t78 = r_i_i_C(3) + pkin(8);
	t57 = qJ(1) + qJ(2);
	t53 = pkin(9) + t57;
	t51 = qJ(4) + t53;
	t47 = sin(t51);
	t48 = cos(t51);
	t59 = cos(qJ(5));
	t70 = qJD(5) * t59;
	t56 = qJD(1) + qJD(2);
	t52 = qJD(4) + t56;
	t58 = sin(qJ(5));
	t74 = t52 * t58;
	t77 = t47 * t74 - t48 * t70;
	t76 = t47 * t70 + t48 * t74;
	t73 = t52 * t59;
	t72 = pkin(1) * qJD(1);
	t71 = qJD(5) * t58;
	t67 = t47 * t71;
	t64 = t47 * t73 + t48 * t71;
	t63 = r_i_i_C(1) * t67 + ((-r_i_i_C(1) * t59 - pkin(4)) * t48 - t78 * t47) * t52 + t76 * r_i_i_C(2);
	t62 = -t64 * r_i_i_C(1) + t77 * r_i_i_C(2) + (-pkin(4) * t47 + t48 * t78) * t52;
	t61 = (-pkin(2) * cos(t57) - pkin(3) * cos(t53)) * t56 + t63;
	t60 = (-pkin(2) * sin(t57) - pkin(3) * sin(t53)) * t56 + t62;
	t1 = [-cos(qJ(1)) * t72 + t61, t61, 0, t63, t77 * r_i_i_C(1) + t64 * r_i_i_C(2); -sin(qJ(1)) * t72 + t60, t60, 0, t62, (-t48 * t73 + t67) * r_i_i_C(2) - t76 * r_i_i_C(1); 0, 0, 0, 0, (-r_i_i_C(1) * t58 - r_i_i_C(2) * t59) * qJD(5);];
	JaD_transl = t1;
end