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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
	% StartTime: 2020-01-03 11:26:17
	% EndTime: 2020-01-03 11:26:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:26:17
	% EndTime: 2020-01-03 11:26:17
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
	% StartTime: 2020-01-03 11:26:17
	% EndTime: 2020-01-03 11:26:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t9 = qJ(1) + pkin(7);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t7 - r_i_i_C(2) * t8) * qJD(1), 0, 0, 0, 0; (cos(qJ(1)) * pkin(1) + r_i_i_C(1) * t8 - r_i_i_C(2) * t7) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:26:17
	% EndTime: 2020-01-03 11:26:17
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (27->11), mult. (32->14), div. (0->0), fcn. (20->6), ass. (0->6)
	t42 = r_i_i_C(3) + qJ(3);
	t41 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(2);
	t38 = qJ(1) + pkin(7);
	t37 = cos(t38);
	t36 = sin(t38);
	t1 = [0, 0, 0, 0, 0; t36 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t42 * t37 - t41 * t36) * qJD(1), 0, qJD(1) * t36, 0, 0; -t37 * qJD(3) + (cos(qJ(1)) * pkin(1) + t42 * t36 + t41 * t37) * qJD(1), 0, -qJD(1) * t37, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:26:17
	% EndTime: 2020-01-03 11:26:17
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (90->26), mult. (128->44), div. (0->0), fcn. (106->8), ass. (0->19)
	t117 = sin(pkin(8));
	t118 = cos(pkin(8));
	t129 = pkin(3) * t118 + pkin(2) + (pkin(6) + r_i_i_C(3)) * t117;
	t119 = sin(qJ(4));
	t127 = t118 * t119;
	t120 = cos(qJ(4));
	t126 = t118 * t120;
	t116 = qJ(1) + pkin(7);
	t114 = sin(t116);
	t115 = cos(t116);
	t124 = t114 * t119 + t115 * t126;
	t123 = t114 * t120 - t115 * t127;
	t122 = -t114 * t126 + t115 * t119;
	t121 = t114 * t127 + t115 * t120;
	t113 = t124 * qJD(1) - t121 * qJD(4);
	t112 = t123 * qJD(1) + t122 * qJD(4);
	t111 = t122 * qJD(1) + t123 * qJD(4);
	t110 = t121 * qJD(1) - t124 * qJD(4);
	t1 = [0, 0, 0, (-r_i_i_C(1) * t120 + r_i_i_C(2) * t119) * t117 * qJD(4), 0; t111 * r_i_i_C(1) + t110 * r_i_i_C(2) + t114 * qJD(3) + (-sin(qJ(1)) * pkin(1) + qJ(3) * t115 - t129 * t114) * qJD(1), 0, qJD(1) * t114, t112 * r_i_i_C(1) - t113 * r_i_i_C(2), 0; t113 * r_i_i_C(1) + t112 * r_i_i_C(2) - t115 * qJD(3) + (cos(qJ(1)) * pkin(1) + qJ(3) * t114 + t129 * t115) * qJD(1), 0, -qJD(1) * t115, -t110 * r_i_i_C(1) + t111 * r_i_i_C(2), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:26:17
	% EndTime: 2020-01-03 11:26:17
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (129->34), mult. (189->54), div. (0->0), fcn. (152->8), ass. (0->25)
	t139 = pkin(4) + r_i_i_C(1);
	t120 = sin(pkin(8));
	t121 = cos(pkin(8));
	t124 = cos(qJ(4));
	t138 = (r_i_i_C(3) + qJ(5) + pkin(6)) * t120 + (t124 * pkin(4) + pkin(3)) * t121 + pkin(2);
	t119 = qJ(1) + pkin(7);
	t117 = sin(t119);
	t118 = cos(t119);
	t123 = sin(qJ(4));
	t136 = t121 * t123;
	t127 = t117 * t136 + t118 * t124;
	t135 = t121 * t124;
	t130 = t117 * t123 + t118 * t135;
	t112 = t127 * qJD(1) - t130 * qJD(4);
	t134 = qJD(1) * t120;
	t133 = qJD(5) * t120;
	t131 = pkin(4) * t123 + qJ(3);
	t129 = t117 * t124 - t118 * t136;
	t128 = -t117 * t135 + t118 * t123;
	t126 = t129 * qJD(4);
	t125 = t127 * qJD(4);
	t114 = t129 * qJD(1) + t128 * qJD(4);
	t115 = t130 * qJD(1) - t125;
	t113 = t128 * qJD(1) + t126;
	t1 = [0, 0, 0, (r_i_i_C(2) * t123 - t139 * t124) * t120 * qJD(4), 0; t118 * t133 + t113 * r_i_i_C(1) + t112 * r_i_i_C(2) + t117 * qJD(3) + pkin(4) * t126 + (-sin(qJ(1)) * pkin(1) + t131 * t118 - t138 * t117) * qJD(1), 0, qJD(1) * t117, -t115 * r_i_i_C(2) + t139 * t114, t118 * t134; t117 * t133 + t115 * r_i_i_C(1) + t114 * r_i_i_C(2) - t118 * qJD(3) - pkin(4) * t125 + (cos(qJ(1)) * pkin(1) + t131 * t117 + t138 * t118) * qJD(1), 0, -qJD(1) * t118, t113 * r_i_i_C(2) - t139 * t112, t117 * t134;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end