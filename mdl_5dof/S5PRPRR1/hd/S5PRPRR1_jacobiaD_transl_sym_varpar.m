% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRPRR1
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:25
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRPRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRPRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:34
	% EndTime: 2019-10-24 10:25:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:34
	% EndTime: 2019-10-24 10:25:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:34
	% EndTime: 2019-10-24 10:25:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (6->3), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->4)
	t32 = pkin(8) + qJ(2);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [0, (-r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(2), 0, 0, 0; 0, (-r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(2), 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:34
	% EndTime: 2019-10-24 10:25:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (24->8), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->6)
	t19 = r_i_i_C(3) + qJ(3);
	t18 = -r_i_i_C(1) * cos(pkin(9)) + r_i_i_C(2) * sin(pkin(9)) - pkin(2);
	t15 = pkin(8) + qJ(2);
	t14 = cos(t15);
	t13 = sin(t15);
	t1 = [0, t14 * qJD(3) + (-t19 * t13 + t18 * t14) * qJD(2), qJD(2) * t14, 0, 0; 0, t13 * qJD(3) + (t18 * t13 + t19 * t14) * qJD(2), qJD(2) * t13, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:34
	% EndTime: 2019-10-24 10:25:34
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (69->20), mult. (70->32), div. (0->0), fcn. (46->5), ass. (0->15)
	t40 = r_i_i_C(3) + pkin(6) + qJ(3);
	t31 = pkin(8) + qJ(2);
	t27 = sin(t31);
	t39 = qJD(2) * t27;
	t29 = cos(t31);
	t38 = qJD(2) * t29;
	t37 = qJD(4) * t27;
	t36 = qJD(4) * t29;
	t30 = pkin(9) + qJ(4);
	t26 = sin(t30);
	t28 = cos(t30);
	t35 = r_i_i_C(1) * t26 + r_i_i_C(2) * t28;
	t34 = -r_i_i_C(1) * t28 + r_i_i_C(2) * t26 - cos(pkin(9)) * pkin(3) - pkin(2);
	t33 = t35 * qJD(4);
	t1 = [0, t29 * qJD(3) + t35 * t37 + (-t40 * t27 + t34 * t29) * qJD(2), t38, (t26 * t36 + t28 * t39) * r_i_i_C(2) + (t26 * t39 - t28 * t36) * r_i_i_C(1), 0; 0, t27 * qJD(3) - t29 * t33 + (t34 * t27 + t40 * t29) * qJD(2), t39, (t26 * t37 - t28 * t38) * r_i_i_C(2) + (-t26 * t38 - t28 * t37) * r_i_i_C(1), 0; 0, 0, 0, -t33, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:34
	% EndTime: 2019-10-24 10:25:34
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (164->31), mult. (120->40), div. (0->0), fcn. (79->7), ass. (0->29)
	t51 = qJD(4) + qJD(5);
	t49 = pkin(9) + qJ(4);
	t47 = qJ(5) + t49;
	t42 = cos(t47);
	t66 = r_i_i_C(2) * t42;
	t41 = sin(t47);
	t68 = r_i_i_C(1) * t41;
	t56 = t66 + t68;
	t54 = t56 * t51;
	t43 = sin(t49);
	t69 = pkin(4) * t43;
	t70 = qJD(4) * t69 + t54;
	t67 = r_i_i_C(2) * t41;
	t65 = r_i_i_C(3) + pkin(7) + pkin(6) + qJ(3);
	t64 = t42 * t51;
	t50 = pkin(8) + qJ(2);
	t44 = sin(t50);
	t63 = qJD(2) * t44;
	t46 = cos(t50);
	t62 = qJD(2) * t46;
	t45 = cos(t49);
	t61 = qJD(4) * t45;
	t60 = r_i_i_C(1) * t64;
	t59 = t51 * t67;
	t57 = qJD(2) * t66;
	t55 = -r_i_i_C(1) * t42 - pkin(4) * t45 - cos(pkin(9)) * pkin(3) - pkin(2) + t67;
	t53 = t44 * t57 + t63 * t68 + (t59 - t60) * t46;
	t36 = t44 * t59;
	t1 = [0, t46 * qJD(3) + t70 * t44 + (-t65 * t44 + t55 * t46) * qJD(2), t62, (t43 * t63 - t46 * t61) * pkin(4) + t53, t53; 0, t44 * qJD(3) - t70 * t46 + (t55 * t44 + t65 * t46) * qJD(2), t63, t36 + (-pkin(4) * t61 - t60) * t44 + (-t56 - t69) * t62, -t46 * t57 + t36 + (-t41 * t62 - t44 * t64) * r_i_i_C(1); 0, 0, 0, -t70, -t54;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end