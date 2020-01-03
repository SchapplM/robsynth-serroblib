% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:52:50
	% EndTime: 2020-01-03 11:52:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:52:50
	% EndTime: 2020-01-03 11:52:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (-r_i_i_C(1) * t5 - r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:52:50
	% EndTime: 2020-01-03 11:52:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t9 = qJ(1) + pkin(9);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t7 - r_i_i_C(2) * t8) * qJD(1), 0, 0, 0, 0; (cos(qJ(1)) * pkin(1) + r_i_i_C(1) * t8 - r_i_i_C(2) * t7) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:52:51
	% EndTime: 2020-01-03 11:52:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (34->9), mult. (24->12), div. (0->0), fcn. (12->6), ass. (0->8)
	t18 = qJ(1) + pkin(9);
	t16 = qJ(3) + t18;
	t14 = sin(t16);
	t15 = cos(t16);
	t17 = qJD(1) + qJD(3);
	t20 = (r_i_i_C(1) * t15 - r_i_i_C(2) * t14) * t17;
	t19 = (-r_i_i_C(1) * t14 - r_i_i_C(2) * t15) * t17;
	t1 = [0, 0, 0, 0, 0; t19 + (-sin(t18) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1), 0, t19, 0, 0; (cos(t18) * pkin(2) + cos(qJ(1)) * pkin(1)) * qJD(1) + t20, 0, t20, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:52:50
	% EndTime: 2020-01-03 11:52:51
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (88->13), mult. (40->15), div. (0->0), fcn. (20->8), ass. (0->13)
	t26 = qJD(1) + qJD(3);
	t32 = pkin(3) * t26;
	t27 = qJ(1) + pkin(9);
	t25 = qJ(3) + t27;
	t23 = qJ(4) + t25;
	t20 = sin(t23);
	t21 = cos(t23);
	t24 = qJD(4) + t26;
	t31 = (r_i_i_C(1) * t21 - r_i_i_C(2) * t20) * t24;
	t30 = cos(t25) * t32 + t31;
	t29 = (-r_i_i_C(1) * t20 - r_i_i_C(2) * t21) * t24;
	t28 = -sin(t25) * t32 + t29;
	t1 = [0, 0, 0, 0, 0; (-sin(t27) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1) + t28, 0, t28, t29, 0; (cos(t27) * pkin(2) + cos(qJ(1)) * pkin(1)) * qJD(1) + t30, 0, t30, t31, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:52:51
	% EndTime: 2020-01-03 11:52:51
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (283->25), mult. (152->38), div. (0->0), fcn. (92->10), ass. (0->26)
	t122 = pkin(8) + r_i_i_C(3);
	t102 = cos(qJ(5));
	t113 = qJD(5) * t102;
	t101 = sin(qJ(5));
	t99 = qJD(1) + qJD(3);
	t97 = qJD(4) + t99;
	t116 = t101 * t97;
	t100 = qJ(1) + pkin(9);
	t98 = qJ(3) + t100;
	t96 = qJ(4) + t98;
	t93 = sin(t96);
	t94 = cos(t96);
	t121 = t94 * t113 - t93 * t116;
	t114 = qJD(5) * t101;
	t115 = t102 * t97;
	t120 = t93 * t114 - t94 * t115;
	t119 = pkin(3) * t99;
	t118 = t93 * t97;
	t117 = t94 * t97;
	t108 = -t93 * t113 - t94 * t116;
	t107 = -t94 * t114 - t93 * t115;
	t106 = pkin(4) * t117 - t120 * r_i_i_C(1) + t108 * r_i_i_C(2) + t122 * t118;
	t105 = cos(t98) * t119 + t106;
	t104 = -pkin(4) * t118 + t107 * r_i_i_C(1) - t121 * r_i_i_C(2) + t122 * t117;
	t103 = -sin(t98) * t119 + t104;
	t1 = [0, 0, 0, 0, (-r_i_i_C(1) * t101 - r_i_i_C(2) * t102) * qJD(5); (-sin(t100) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1) + t103, 0, t103, t104, t108 * r_i_i_C(1) + t120 * r_i_i_C(2); (cos(t100) * pkin(2) + cos(qJ(1)) * pkin(1)) * qJD(1) + t105, 0, t105, t106, t121 * r_i_i_C(1) + t107 * r_i_i_C(2);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end