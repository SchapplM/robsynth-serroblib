% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:20
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRPR10_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR10_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t81 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t81, 0, 0, 0, 0; 0, sin(qJ(1)) * t81, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t108 = sin(qJ(2));
	t109 = sin(qJ(1));
	t116 = t108 * t109;
	t111 = cos(qJ(1));
	t115 = t108 * t111;
	t110 = cos(qJ(2));
	t114 = t109 * t110;
	t113 = t110 * t111;
	t106 = sin(pkin(6));
	t112 = qJD(1) * t106;
	t107 = cos(pkin(6));
	t1 = [0, t111 * t112, 0, (-t107 * t116 + t113) * qJD(2) + (t107 * t113 - t116) * qJD(1), 0, 0; 0, t109 * t112, 0, (t107 * t115 + t114) * qJD(2) + (t107 * t114 + t115) * qJD(1), 0, 0; 0, 0, 0, t106 * qJD(2) * t108, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:53
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t124 = sin(qJ(2));
	t125 = sin(qJ(1));
	t132 = t124 * t125;
	t127 = cos(qJ(1));
	t131 = t124 * t127;
	t126 = cos(qJ(2));
	t130 = t125 * t126;
	t129 = t126 * t127;
	t122 = sin(pkin(6));
	t128 = qJD(1) * t122;
	t123 = cos(pkin(6));
	t1 = [0, t127 * t128, 0, (-t123 * t132 + t129) * qJD(2) + (t123 * t129 - t132) * qJD(1), 0, 0; 0, t125 * t128, 0, (t123 * t131 + t130) * qJD(2) + (t123 * t130 + t131) * qJD(1), 0, 0; 0, 0, 0, t122 * qJD(2) * t124, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:53
	% EndTime: 2019-10-10 10:20:53
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (33->17), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->23)
	t163 = sin(pkin(6));
	t166 = sin(qJ(1));
	t181 = t163 * t166;
	t168 = cos(qJ(1));
	t180 = t163 * t168;
	t165 = sin(qJ(2));
	t179 = t165 * t166;
	t178 = t165 * t168;
	t167 = cos(qJ(2));
	t177 = t166 * t167;
	t176 = t168 * t167;
	t175 = qJD(1) * t163;
	t162 = pkin(11) + qJ(4);
	t161 = cos(t162);
	t174 = qJD(2) * t161;
	t173 = qJD(2) * t163;
	t164 = cos(pkin(6));
	t172 = t164 * t176 - t179;
	t171 = t164 * t177 + t178;
	t170 = t164 * t178 + t177;
	t169 = -t164 * t179 + t176;
	t160 = sin(t162);
	t1 = [0, t168 * t175, 0, t172 * qJD(1) + t169 * qJD(2), 0, (-t169 * t160 + t161 * t181) * qJD(4) - t171 * t174 + (t160 * t180 - t170 * t161) * qJD(1); 0, t166 * t175, 0, t171 * qJD(1) + t170 * qJD(2), 0, (-t170 * t160 - t161 * t180) * qJD(4) + t172 * t174 + (t160 * t181 + t169 * t161) * qJD(1); 0, 0, 0, t165 * t173, 0, t167 * t161 * t173 + (-t160 * t163 * t165 + t161 * t164) * qJD(4);];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end