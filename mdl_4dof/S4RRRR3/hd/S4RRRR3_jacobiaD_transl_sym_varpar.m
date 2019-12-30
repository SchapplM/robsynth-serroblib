% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4RRRR3
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4RRRR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RRRR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRRR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:07:06
	% EndTime: 2019-12-29 15:07:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:07:06
	% EndTime: 2019-12-29 15:07:06
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:07:06
	% EndTime: 2019-12-29 15:07:06
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(5) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t18 * t28 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0; -t20 * t21 + (t22 * t18 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0; 0, -t21, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:07:06
	% EndTime: 2019-12-29 15:07:06
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (81->26), mult. (114->37), div. (0->0), fcn. (73->6), ass. (0->27)
	t38 = qJD(2) + qJD(3);
	t39 = qJ(2) + qJ(3);
	t37 = cos(t39);
	t59 = r_i_i_C(2) * t37;
	t36 = sin(t39);
	t61 = r_i_i_C(1) * t36;
	t49 = t59 + t61;
	t47 = t49 * t38;
	t40 = sin(qJ(2));
	t62 = pkin(2) * t40;
	t63 = qJD(2) * t62 + t47;
	t60 = r_i_i_C(2) * t36;
	t58 = r_i_i_C(3) + pkin(6) + pkin(5);
	t57 = t37 * t38;
	t41 = sin(qJ(1));
	t56 = qJD(1) * t41;
	t43 = cos(qJ(1));
	t55 = qJD(1) * t43;
	t42 = cos(qJ(2));
	t54 = qJD(2) * t42;
	t53 = r_i_i_C(1) * t57;
	t52 = t38 * t60;
	t51 = qJD(1) * t59;
	t48 = -t42 * pkin(2) - r_i_i_C(1) * t37 - pkin(1) + t60;
	t46 = t41 * t51 + t56 * t61 + (t52 - t53) * t43;
	t31 = t41 * t52;
	t1 = [t63 * t41 + (-t41 * t58 + t43 * t48) * qJD(1), (t40 * t56 - t43 * t54) * pkin(2) + t46, t46, 0; -t63 * t43 + (t41 * t48 + t43 * t58) * qJD(1), t31 + (-pkin(2) * t54 - t53) * t41 + (-t49 - t62) * t55, -t43 * t51 + t31 + (-t36 * t55 - t41 * t57) * r_i_i_C(1), 0; 0, -t63, -t47, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:07:06
	% EndTime: 2019-12-29 15:07:06
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (202->36), mult. (170->47), div. (0->0), fcn. (109->8), ass. (0->37)
	t54 = sin(qJ(2));
	t51 = qJD(2) + qJD(3);
	t47 = qJD(4) + t51;
	t53 = qJ(2) + qJ(3);
	t50 = qJ(4) + t53;
	t46 = cos(t50);
	t77 = r_i_i_C(2) * t46;
	t45 = sin(t50);
	t79 = r_i_i_C(1) * t45;
	t62 = t77 + t79;
	t60 = t62 * t47;
	t48 = sin(t53);
	t80 = pkin(3) * t48;
	t58 = -t51 * t80 - t60;
	t73 = pkin(2) * qJD(2);
	t82 = -t54 * t73 + t58;
	t75 = t46 * t47;
	t68 = r_i_i_C(1) * t75;
	t49 = cos(t53);
	t74 = t49 * t51;
	t81 = -pkin(3) * t74 - t68;
	t78 = r_i_i_C(2) * t45;
	t76 = r_i_i_C(3) + pkin(7) + pkin(6) + pkin(5);
	t55 = sin(qJ(1));
	t72 = qJD(1) * t55;
	t57 = cos(qJ(1));
	t71 = qJD(1) * t57;
	t67 = t47 * t78;
	t65 = qJD(1) * t77;
	t66 = t55 * t65 + t57 * t67 + t72 * t79;
	t56 = cos(qJ(2));
	t63 = -t56 * t73 + t81;
	t61 = -t56 * pkin(2) - pkin(3) * t49 - r_i_i_C(1) * t46 - pkin(1) + t78;
	t59 = -t57 * t68 + t66;
	t44 = -t54 * pkin(2) - t80;
	t39 = t55 * t67;
	t1 = [-t82 * t55 + (-t76 * t55 + t61 * t57) * qJD(1), -t44 * t72 + t63 * t57 + t66, (t48 * t72 - t57 * t74) * pkin(3) + t59, t59; t82 * t57 + (t61 * t55 + t76 * t57) * qJD(1), t39 + t63 * t55 + (t44 - t62) * t71, t39 + t81 * t55 + (-t62 - t80) * t71, -t57 * t65 + t39 + (-t45 * t71 - t55 * t75) * r_i_i_C(1); 0, t82, t58, -t60;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end