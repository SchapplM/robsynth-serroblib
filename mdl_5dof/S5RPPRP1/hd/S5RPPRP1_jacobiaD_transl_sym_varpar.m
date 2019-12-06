% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPRP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPRP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:36:51
	% EndTime: 2019-12-05 17:36:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:36:51
	% EndTime: 2019-12-05 17:36:51
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
	% StartTime: 2019-12-05 17:36:51
	% EndTime: 2019-12-05 17:36:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t9 = qJ(1) + pkin(7);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, 0, 0, 0; (sin(qJ(1)) * pkin(1) + r_i_i_C(1) * t7 + r_i_i_C(2) * t8) * qJD(1), 0, 0, 0, 0; (-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t8 + r_i_i_C(2) * t7) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:36:51
	% EndTime: 2019-12-05 17:36:51
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (27->11), mult. (32->14), div. (0->0), fcn. (20->6), ass. (0->6)
	t42 = -r_i_i_C(3) - qJ(3);
	t41 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(2);
	t38 = qJ(1) + pkin(7);
	t37 = cos(t38);
	t36 = sin(t38);
	t1 = [0, 0, 0, 0, 0; -t36 * qJD(3) + (sin(qJ(1)) * pkin(1) + t42 * t37 + t41 * t36) * qJD(1), 0, -qJD(1) * t36, 0, 0; t37 * qJD(3) + (-cos(qJ(1)) * pkin(1) + t42 * t36 - t41 * t37) * qJD(1), 0, qJD(1) * t37, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:36:51
	% EndTime: 2019-12-05 17:36:51
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (90->26), mult. (128->44), div. (0->0), fcn. (106->8), ass. (0->19)
	t121 = sin(pkin(8));
	t122 = cos(pkin(8));
	t133 = pkin(3) * t122 + pkin(2) + (pkin(6) + r_i_i_C(3)) * t121;
	t123 = sin(qJ(4));
	t131 = t122 * t123;
	t124 = cos(qJ(4));
	t130 = t122 * t124;
	t120 = qJ(1) + pkin(7);
	t118 = sin(t120);
	t119 = cos(t120);
	t128 = t118 * t123 + t119 * t130;
	t127 = -t118 * t124 + t119 * t131;
	t126 = t118 * t130 - t119 * t123;
	t125 = t118 * t131 + t119 * t124;
	t117 = t128 * qJD(1) - t125 * qJD(4);
	t116 = t127 * qJD(1) + t126 * qJD(4);
	t115 = t126 * qJD(1) + t127 * qJD(4);
	t114 = t125 * qJD(1) - t128 * qJD(4);
	t1 = [0, 0, 0, (-r_i_i_C(1) * t124 + r_i_i_C(2) * t123) * t121 * qJD(4), 0; t115 * r_i_i_C(1) - t114 * r_i_i_C(2) - t118 * qJD(3) + (sin(qJ(1)) * pkin(1) - qJ(3) * t119 + t133 * t118) * qJD(1), 0, -qJD(1) * t118, t116 * r_i_i_C(1) + t117 * r_i_i_C(2), 0; -t117 * r_i_i_C(1) + t116 * r_i_i_C(2) + t119 * qJD(3) + (-cos(qJ(1)) * pkin(1) - qJ(3) * t118 - t133 * t119) * qJD(1), 0, qJD(1) * t119, t114 * r_i_i_C(1) + t115 * r_i_i_C(2), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:36:51
	% EndTime: 2019-12-05 17:36:51
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (131->36), mult. (189->54), div. (0->0), fcn. (152->8), ass. (0->25)
	t143 = pkin(4) + r_i_i_C(1);
	t123 = qJ(1) + pkin(7);
	t121 = sin(t123);
	t122 = cos(t123);
	t128 = cos(qJ(4));
	t125 = cos(pkin(8));
	t127 = sin(qJ(4));
	t139 = t125 * t127;
	t130 = t121 * t139 + t122 * t128;
	t142 = t130 * qJD(4);
	t124 = sin(pkin(8));
	t141 = (r_i_i_C(3) + qJ(5) + pkin(6)) * t124 + (t128 * pkin(4) + pkin(3)) * t125 + pkin(2);
	t138 = t125 * t128;
	t137 = qJD(1) * t124;
	t136 = qJD(5) * t124;
	t134 = -pkin(4) * t127 - qJ(3);
	t133 = t121 * t127 + t122 * t138;
	t132 = -t121 * t128 + t122 * t139;
	t131 = t121 * t138 - t122 * t127;
	t129 = t132 * qJD(4);
	t116 = t130 * qJD(1) - t133 * qJD(4);
	t118 = t132 * qJD(1) + t131 * qJD(4);
	t119 = t133 * qJD(1) - t142;
	t117 = t131 * qJD(1) + t129;
	t1 = [0, 0, 0, (r_i_i_C(2) * t127 - t143 * t128) * t124 * qJD(4), 0; -t122 * t136 + t117 * r_i_i_C(1) - t116 * r_i_i_C(2) - t121 * qJD(3) + pkin(4) * t129 + (sin(qJ(1)) * pkin(1) + t134 * t122 + t141 * t121) * qJD(1), 0, -qJD(1) * t121, t119 * r_i_i_C(2) + t143 * t118, -t122 * t137; -t121 * t136 - t119 * r_i_i_C(1) + t118 * r_i_i_C(2) + t122 * qJD(3) + pkin(4) * t142 + (-cos(qJ(1)) * pkin(1) + t134 * t121 - t141 * t122) * qJD(1), 0, qJD(1) * t122, t117 * r_i_i_C(2) + t143 * t116, -t121 * t137;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end