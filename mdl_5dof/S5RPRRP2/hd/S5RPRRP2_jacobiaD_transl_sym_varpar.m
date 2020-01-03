% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRP2
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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRP2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRP2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:45:53
	% EndTime: 2020-01-03 11:45:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:45:53
	% EndTime: 2020-01-03 11:45:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (-r_i_i_C(1) * t5 - r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:45:53
	% EndTime: 2020-01-03 11:45:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t9 = qJ(1) + pkin(8);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t7 - r_i_i_C(2) * t8) * qJD(1), 0, 0, 0, 0; (cos(qJ(1)) * pkin(1) + r_i_i_C(1) * t8 - r_i_i_C(2) * t7) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:45:53
	% EndTime: 2020-01-03 11:45:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (34->9), mult. (24->12), div. (0->0), fcn. (12->6), ass. (0->8)
	t18 = qJ(1) + pkin(8);
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
	% StartTime: 2020-01-03 11:45:53
	% EndTime: 2020-01-03 11:45:53
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (133->21), mult. (108->35), div. (0->0), fcn. (66->8), ass. (0->21)
	t110 = pkin(7) + r_i_i_C(3);
	t93 = cos(qJ(4));
	t102 = qJD(4) * t93;
	t90 = qJD(1) + qJD(3);
	t92 = sin(qJ(4));
	t105 = t90 * t92;
	t91 = qJ(1) + pkin(8);
	t89 = qJ(3) + t91;
	t87 = sin(t89);
	t88 = cos(t89);
	t109 = t88 * t102 - t87 * t105;
	t103 = qJD(4) * t92;
	t104 = t90 * t93;
	t108 = t87 * t103 - t88 * t104;
	t107 = t87 * t90;
	t106 = t88 * t90;
	t97 = -t88 * t103 - t87 * t104;
	t96 = -t87 * t102 - t88 * t105;
	t95 = pkin(3) * t106 - t108 * r_i_i_C(1) + t96 * r_i_i_C(2) + t107 * t110;
	t94 = -pkin(3) * t107 + t97 * r_i_i_C(1) - t109 * r_i_i_C(2) + t106 * t110;
	t1 = [0, 0, 0, (-r_i_i_C(1) * t92 - r_i_i_C(2) * t93) * qJD(4), 0; (-sin(qJ(1)) * pkin(1) - sin(t91) * pkin(2)) * qJD(1) + t94, 0, t94, t96 * r_i_i_C(1) + r_i_i_C(2) * t108, 0; (cos(t91) * pkin(2) + cos(qJ(1)) * pkin(1)) * qJD(1) + t95, 0, t95, r_i_i_C(1) * t109 + t97 * r_i_i_C(2), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:45:53
	% EndTime: 2020-01-03 11:45:53
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (187->23), mult. (144->29), div. (0->0), fcn. (93->8), ass. (0->22)
	t101 = cos(qJ(4));
	t113 = pkin(4) + r_i_i_C(1);
	t118 = t113 * t101;
	t117 = pkin(3) + t118;
	t100 = sin(qJ(4));
	t111 = r_i_i_C(2) * t100;
	t116 = qJD(4) * (-t111 + t118);
	t97 = qJD(1) + qJD(3);
	t115 = t97 * t111 + qJD(5);
	t98 = qJ(1) + pkin(8);
	t96 = qJ(3) + t98;
	t93 = sin(t96);
	t110 = t97 * t93;
	t94 = cos(t96);
	t109 = t97 * t94;
	t107 = -r_i_i_C(2) * t101 - t113 * t100;
	t106 = t107 * t97;
	t105 = t107 * qJD(4);
	t104 = -t97 * (-qJ(5) - pkin(7)) + t105;
	t103 = r_i_i_C(3) * t109 + t104 * t94 - t117 * t110 + t115 * t93;
	t102 = r_i_i_C(3) * t110 + t104 * t93 + t117 * t109 - t115 * t94;
	t1 = [0, 0, 0, t105, 0; (-sin(t98) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1) + t103, 0, t103, t94 * t106 - t93 * t116, t110; (cos(t98) * pkin(2) + cos(qJ(1)) * pkin(1)) * qJD(1) + t102, 0, t102, t93 * t106 + t94 * t116, -t109;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end