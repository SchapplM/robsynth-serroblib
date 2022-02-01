% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
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
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(8);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (43->18), mult. (68->31), div. (0->0), fcn. (42->6), ass. (0->13)
	t24 = sin(qJ(3));
	t25 = cos(qJ(3));
	t26 = (r_i_i_C(1) * t24 + r_i_i_C(2) * t25) * qJD(3);
	t33 = pkin(6) + r_i_i_C(3);
	t32 = qJD(1) * t24;
	t31 = qJD(1) * t25;
	t30 = qJD(3) * t24;
	t29 = qJD(3) * t25;
	t27 = -r_i_i_C(1) * t25 + r_i_i_C(2) * t24 - pkin(2);
	t23 = qJ(1) + pkin(8);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [t21 * t26 + (-cos(qJ(1)) * pkin(1) - t33 * t21 + t27 * t22) * qJD(1), 0, (t21 * t31 + t22 * t30) * r_i_i_C(2) + (t21 * t32 - t22 * t29) * r_i_i_C(1), 0, 0; -t22 * t26 + (-sin(qJ(1)) * pkin(1) + t33 * t22 + t27 * t21) * qJD(1), 0, (t21 * t30 - t22 * t31) * r_i_i_C(2) + (-t21 * t29 - t22 * t32) * r_i_i_C(1), 0, 0; 0, 0, -t26, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (119->29), mult. (118->40), div. (0->0), fcn. (75->8), ass. (0->26)
	t44 = qJD(3) + qJD(4);
	t46 = qJ(3) + qJ(4);
	t42 = sin(t46);
	t43 = cos(t46);
	t63 = r_i_i_C(2) * t43;
	t54 = r_i_i_C(1) * t42 + t63;
	t52 = t54 * t44;
	t47 = sin(qJ(3));
	t65 = pkin(3) * t47;
	t66 = qJD(3) * t65 + t52;
	t64 = r_i_i_C(2) * t42;
	t62 = r_i_i_C(3) + pkin(7) + pkin(6);
	t61 = t43 * t44;
	t60 = qJD(1) * t42;
	t48 = cos(qJ(3));
	t59 = qJD(3) * t48;
	t58 = r_i_i_C(1) * t61;
	t57 = t44 * t64;
	t56 = qJD(1) * t63;
	t53 = -t48 * pkin(3) - r_i_i_C(1) * t43 - pkin(2) + t64;
	t45 = qJ(1) + pkin(8);
	t40 = sin(t45);
	t41 = cos(t45);
	t51 = (t57 - t58) * t41 + (t60 * r_i_i_C(1) + t56) * t40;
	t35 = t40 * t57;
	t1 = [t66 * t40 + (-cos(qJ(1)) * pkin(1) - t62 * t40 + t53 * t41) * qJD(1), 0, (qJD(1) * t40 * t47 - t41 * t59) * pkin(3) + t51, t51, 0; -t66 * t41 + (-sin(qJ(1)) * pkin(1) + t62 * t41 + t53 * t40) * qJD(1), 0, t35 + (-pkin(3) * t59 - t58) * t40 + (-t54 - t65) * t41 * qJD(1), -t41 * t56 + t35 + (-t40 * t61 - t41 * t60) * r_i_i_C(1), 0; 0, 0, -t66, -t52, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (167->35), mult. (150->46), div. (0->0), fcn. (97->8), ass. (0->32)
	t53 = qJ(3) + qJ(4);
	t49 = cos(t53);
	t51 = qJD(3) + qJD(4);
	t69 = t49 * t51;
	t73 = -pkin(4) - r_i_i_C(1);
	t75 = t73 * t69;
	t48 = sin(t53);
	t72 = r_i_i_C(2) * t49;
	t59 = r_i_i_C(1) * t48 + t72;
	t54 = sin(qJ(3));
	t68 = pkin(3) * qJD(3);
	t62 = t54 * t68;
	t70 = t48 * t51;
	t74 = pkin(4) * t70 + t59 * t51 + t62;
	t71 = r_i_i_C(3) + qJ(5) + pkin(7) + pkin(6);
	t52 = qJ(1) + pkin(8);
	t46 = sin(t52);
	t67 = qJD(1) * t46;
	t47 = cos(t52);
	t66 = qJD(1) * t47;
	t65 = r_i_i_C(2) * t70;
	t64 = t47 * t69;
	t61 = t48 * t67;
	t63 = r_i_i_C(1) * t61 + t47 * t65 + t67 * t72;
	t55 = cos(qJ(3));
	t60 = -t55 * t68 + t75;
	t58 = -pkin(3) * t55 + r_i_i_C(2) * t48 + t73 * t49 - pkin(2);
	t57 = t73 * t48 - t72;
	t56 = t57 * t51;
	t45 = -pkin(3) * t54 - pkin(4) * t48;
	t40 = t46 * t65;
	t1 = [t47 * qJD(5) + t74 * t46 + (-cos(qJ(1)) * pkin(1) - t71 * t46 + t58 * t47) * qJD(1), 0, -t45 * t67 + t60 * t47 + t63, -r_i_i_C(1) * t64 + (t61 - t64) * pkin(4) + t63, t66; t46 * qJD(5) - t74 * t47 + (-sin(qJ(1)) * pkin(1) + t71 * t47 + t58 * t46) * qJD(1), 0, t40 + t60 * t46 + (t45 - t59) * t66, t46 * t75 + t57 * t66 + t40, t67; 0, 0, t56 - t62, t56, 0;];
	JaD_transl = t1;
end