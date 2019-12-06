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
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
	% StartTime: 2019-12-05 18:02:09
	% EndTime: 2019-12-05 18:02:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:02:09
	% EndTime: 2019-12-05 18:02:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (r_i_i_C(1) * t5 + r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t6 + r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:02:09
	% EndTime: 2019-12-05 18:02:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t9 = qJ(1) + pkin(8);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, 0, 0, 0; (sin(qJ(1)) * pkin(1) + r_i_i_C(1) * t7 + r_i_i_C(2) * t8) * qJD(1), 0, 0, 0, 0; (-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t8 + r_i_i_C(2) * t7) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:02:09
	% EndTime: 2019-12-05 18:02:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (34->9), mult. (24->12), div. (0->0), fcn. (12->6), ass. (0->8)
	t18 = qJ(1) + pkin(8);
	t16 = qJ(3) + t18;
	t17 = qJD(1) + qJD(3);
	t22 = sin(t16) * t17;
	t21 = cos(t16) * t17;
	t20 = r_i_i_C(1) * t22 + r_i_i_C(2) * t21;
	t19 = -r_i_i_C(1) * t21 + r_i_i_C(2) * t22;
	t1 = [0, 0, 0, 0, 0; (sin(t18) * pkin(2) + sin(qJ(1)) * pkin(1)) * qJD(1) + t20, 0, t20, 0, 0; (-cos(t18) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t19, 0, t19, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:02:09
	% EndTime: 2019-12-05 18:02:09
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (133->23), mult. (108->37), div. (0->0), fcn. (66->8), ass. (0->19)
	t92 = cos(qJ(4));
	t101 = qJD(4) * t92;
	t89 = qJD(1) + qJD(3);
	t91 = sin(qJ(4));
	t104 = t89 * t91;
	t90 = qJ(1) + pkin(8);
	t88 = qJ(3) + t90;
	t86 = sin(t88);
	t87 = cos(t88);
	t107 = t86 * t101 + t87 * t104;
	t102 = qJD(4) * t91;
	t103 = t89 * t92;
	t106 = t87 * t102 + t86 * t103;
	t105 = -pkin(7) - r_i_i_C(3);
	t98 = t86 * t102;
	t95 = t87 * t101;
	t94 = r_i_i_C(2) * t95 + (t105 * t87 + (-r_i_i_C(2) * t91 + pkin(3)) * t86) * t89 + t106 * r_i_i_C(1);
	t93 = r_i_i_C(1) * t98 + ((-r_i_i_C(1) * t92 - pkin(3)) * t87 + t105 * t86) * t89 + t107 * r_i_i_C(2);
	t1 = [0, 0, 0, (-r_i_i_C(1) * t91 - r_i_i_C(2) * t92) * qJD(4), 0; (sin(t90) * pkin(2) + sin(qJ(1)) * pkin(1)) * qJD(1) + t94, 0, t94, (t87 * t103 - t98) * r_i_i_C(2) + t107 * r_i_i_C(1), 0; (-cos(t90) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t93, 0, t93, t106 * r_i_i_C(2) + (t86 * t104 - t95) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:02:09
	% EndTime: 2019-12-05 18:02:09
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (187->24), mult. (144->28), div. (0->0), fcn. (93->8), ass. (0->21)
	t105 = sin(qJ(4));
	t106 = cos(qJ(4));
	t120 = pkin(4) + r_i_i_C(1);
	t122 = r_i_i_C(2) * t106 + t120 * t105;
	t130 = t122 * qJD(4);
	t129 = t120 * t106;
	t102 = qJD(1) + qJD(3);
	t127 = (-r_i_i_C(3) - qJ(5) - pkin(7)) * t102 + t130;
	t126 = pkin(3) + t129;
	t118 = r_i_i_C(2) * t105;
	t125 = qJD(4) * (-t118 + t129);
	t103 = qJ(1) + pkin(8);
	t101 = qJ(3) + t103;
	t98 = sin(t101);
	t116 = t102 * t98;
	t99 = cos(t101);
	t115 = t102 * t99;
	t109 = t102 * t122;
	t108 = t115 * t118 + (-t126 * t102 + qJD(5)) * t99 + t127 * t98;
	t107 = (-t118 * t102 - qJD(5)) * t98 + t126 * t116 + t127 * t99;
	t1 = [0, 0, 0, -t130, 0; (sin(t103) * pkin(2) + sin(qJ(1)) * pkin(1)) * qJD(1) + t107, 0, t107, t99 * t109 + t98 * t125, -t116; (-cos(t103) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t108, 0, t108, t98 * t109 - t99 * t125, t115;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end