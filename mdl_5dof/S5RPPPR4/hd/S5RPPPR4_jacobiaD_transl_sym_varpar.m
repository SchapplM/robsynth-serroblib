% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPPR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPPR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:48:47
	% EndTime: 2019-12-29 15:48:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:48:47
	% EndTime: 2019-12-29 15:48:47
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:48:47
	% EndTime: 2019-12-29 15:48:47
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(7);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:48:43
	% EndTime: 2019-12-29 15:48:43
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (22->9), mult. (24->12), div. (0->0), fcn. (14->4), ass. (0->6)
	t13 = -pkin(2) + r_i_i_C(2);
	t12 = r_i_i_C(3) + qJ(3);
	t11 = qJ(1) + pkin(7);
	t10 = cos(t11);
	t9 = sin(t11);
	t1 = [t10 * qJD(3) + (-cos(qJ(1)) * pkin(1) - t12 * t9 + t13 * t10) * qJD(1), 0, qJD(1) * t10, 0, 0; t9 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t13 * t9 + t12 * t10) * qJD(1), 0, qJD(1) * t9, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:48:52
	% EndTime: 2019-12-29 15:48:52
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (37->13), mult. (40->16), div. (0->0), fcn. (26->6), ass. (0->8)
	t16 = qJ(1) + pkin(7);
	t14 = sin(t16);
	t21 = qJD(1) * t14;
	t20 = -pkin(2) - r_i_i_C(3) - qJ(4);
	t19 = r_i_i_C(1) * sin(pkin(8)) + r_i_i_C(2) * cos(pkin(8)) + qJ(3);
	t15 = cos(t16);
	t13 = qJD(1) * t15;
	t1 = [t15 * qJD(3) - t14 * qJD(4) + (-cos(qJ(1)) * pkin(1) + t20 * t15 - t19 * t14) * qJD(1), 0, t13, -t21, 0; t14 * qJD(3) + t15 * qJD(4) + (-sin(qJ(1)) * pkin(1) + t19 * t15 + t20 * t14) * qJD(1), 0, t21, t13, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:48:42
	% EndTime: 2019-12-29 15:48:42
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (84->24), mult. (86->34), div. (0->0), fcn. (56->7), ass. (0->15)
	t30 = pkin(8) + qJ(5);
	t26 = sin(t30);
	t28 = cos(t30);
	t41 = (r_i_i_C(1) * t28 - r_i_i_C(2) * t26) * qJD(5);
	t31 = qJ(1) + pkin(7);
	t27 = sin(t31);
	t40 = qJD(1) * t27;
	t29 = cos(t31);
	t25 = qJD(1) * t29;
	t39 = qJD(5) * t27;
	t38 = qJD(5) * t29;
	t37 = -pkin(2) - r_i_i_C(3) - pkin(6) - qJ(4);
	t35 = pkin(4) * sin(pkin(8)) + r_i_i_C(1) * t26 + r_i_i_C(2) * t28 + qJ(3);
	t34 = qJD(3) + t41;
	t1 = [-t27 * qJD(4) + t34 * t29 + (-cos(qJ(1)) * pkin(1) + t37 * t29 - t35 * t27) * qJD(1), 0, t25, -t40, (-t26 * t25 - t28 * t39) * r_i_i_C(2) + (t28 * t25 - t26 * t39) * r_i_i_C(1); t29 * qJD(4) + t34 * t27 + (-sin(qJ(1)) * pkin(1) + t37 * t27 + t35 * t29) * qJD(1), 0, t40, t25, (-t26 * t40 + t28 * t38) * r_i_i_C(2) + (t26 * t38 + t28 * t40) * r_i_i_C(1); 0, 0, 0, 0, -t41;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end