% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:03
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:03:38
	% EndTime: 2020-11-04 22:03:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:03:38
	% EndTime: 2020-11-04 22:03:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t81 = cos(qJ(1));
	t80 = sin(qJ(1));
	t1 = [t81, -t80, 0, 0; t80, t81, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:03:38
	% EndTime: 2020-11-04 22:03:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t82 = sin(pkin(6));
	t85 = sin(qJ(1));
	t93 = t85 * t82;
	t84 = sin(qJ(2));
	t92 = t85 * t84;
	t86 = cos(qJ(2));
	t91 = t85 * t86;
	t87 = cos(qJ(1));
	t90 = t87 * t82;
	t89 = t87 * t84;
	t88 = t87 * t86;
	t83 = cos(pkin(6));
	t1 = [-t83 * t92 + t88, -t83 * t91 - t89, t93, t87 * pkin(1) + pkin(8) * t93 + 0; t83 * t89 + t91, t83 * t88 - t92, -t90, t85 * pkin(1) - pkin(8) * t90 + 0; t82 * t84, t82 * t86, t83, t83 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:03:38
	% EndTime: 2020-11-04 22:03:38
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (22->17), mult. (39->24), div. (0->0), fcn. (52->6), ass. (0->14)
	t98 = sin(qJ(2));
	t99 = sin(qJ(1));
	t106 = t99 * t98;
	t101 = cos(qJ(1));
	t105 = t101 * t98;
	t100 = cos(qJ(2));
	t104 = t99 * t100;
	t103 = t101 * t100;
	t102 = pkin(2) * t98 - qJ(3) * t100;
	t97 = cos(pkin(6));
	t96 = sin(pkin(6));
	t95 = t100 * pkin(2) + t98 * qJ(3) + pkin(1);
	t94 = t96 * pkin(8) - t102 * t97;
	t1 = [-t97 * t106 + t103, t99 * t96, t97 * t104 + t105, t95 * t101 + t94 * t99 + 0; t97 * t105 + t104, -t101 * t96, -t97 * t103 + t106, -t94 * t101 + t95 * t99 + 0; t96 * t98, t97, -t96 * t100, t97 * pkin(8) + t102 * t96 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:03:38
	% EndTime: 2020-11-04 22:03:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->20), mult. (39->24), div. (0->0), fcn. (52->6), ass. (0->16)
	t112 = sin(qJ(2));
	t113 = sin(qJ(1));
	t121 = t113 * t112;
	t114 = cos(qJ(2));
	t120 = t113 * t114;
	t115 = cos(qJ(1));
	t119 = t115 * t112;
	t118 = t115 * t114;
	t116 = pkin(2) + pkin(3);
	t117 = qJ(3) * t114 - t112 * t116;
	t111 = qJ(4) - pkin(8);
	t110 = cos(pkin(6));
	t109 = sin(pkin(6));
	t108 = t112 * qJ(3) + t116 * t114 + pkin(1);
	t107 = -t109 * t111 + t117 * t110;
	t1 = [-t110 * t121 + t118, t110 * t120 + t119, -t113 * t109, t107 * t113 + t108 * t115 + 0; t110 * t119 + t120, -t110 * t118 + t121, t115 * t109, -t107 * t115 + t108 * t113 + 0; t109 * t112, -t109 * t114, -t110, -t117 * t109 - t111 * t110 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:03:39
	% EndTime: 2020-11-04 22:03:39
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (46->25), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->23)
	t126 = sin(pkin(6));
	t133 = cos(qJ(5));
	t143 = t126 * t133;
	t130 = sin(qJ(5));
	t142 = t130 * t126;
	t131 = sin(qJ(2));
	t141 = t131 * t133;
	t132 = sin(qJ(1));
	t140 = t132 * t131;
	t134 = cos(qJ(2));
	t139 = t132 * t134;
	t135 = cos(qJ(1));
	t138 = t135 * t131;
	t137 = t135 * t134;
	t125 = pkin(2) + pkin(3) + pkin(4);
	t129 = qJ(3) - pkin(9);
	t136 = t125 * t131 - t129 * t134;
	t128 = qJ(4) - pkin(8);
	t127 = cos(pkin(6));
	t124 = t125 * t134 + t129 * t131 + pkin(1);
	t123 = t127 * t138 + t139;
	t122 = -t126 * t128 - t136 * t127;
	t1 = [(-t127 * t141 - t142) * t132 + t133 * t137, (t127 * t140 - t137) * t130 - t132 * t143, -t127 * t139 - t138, t122 * t132 + t124 * t135 + 0; t123 * t133 + t135 * t142, -t123 * t130 + t135 * t143, t127 * t137 - t140, -t122 * t135 + t124 * t132 + 0; t126 * t141 - t127 * t130, -t127 * t133 - t131 * t142, t126 * t134, t136 * t126 - t128 * t127 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:03:39
	% EndTime: 2020-11-04 22:03:39
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (75->38), mult. (112->60), div. (0->0), fcn. (146->10), ass. (0->30)
	t150 = sin(pkin(6));
	t158 = cos(qJ(5));
	t172 = t150 * t158;
	t153 = sin(qJ(6));
	t159 = cos(qJ(2));
	t171 = t153 * t159;
	t154 = sin(qJ(5));
	t170 = t154 * t150;
	t155 = sin(qJ(2));
	t169 = t155 * t158;
	t156 = sin(qJ(1));
	t168 = t156 * t159;
	t157 = cos(qJ(6));
	t167 = t157 * t159;
	t166 = t158 * t159;
	t160 = cos(qJ(1));
	t165 = t160 * t155;
	t164 = t160 * t159;
	t149 = pkin(5) * t158 + pkin(10) * t154 + pkin(2) + pkin(3) + pkin(4);
	t152 = qJ(3) - pkin(9);
	t163 = t149 * t155 - t152 * t159;
	t151 = cos(pkin(6));
	t147 = t151 * t169 + t170;
	t162 = -t147 * t153 + t151 * t167;
	t161 = pkin(5) * t154 - pkin(10) * t158 - pkin(8) + qJ(4);
	t148 = t153 * t166 + t155 * t157;
	t146 = t150 * t169 - t151 * t154;
	t145 = t149 * t159 + t152 * t155 + pkin(1);
	t144 = -t150 * t161 - t151 * t163;
	t1 = [(-t147 * t156 + t158 * t164) * t157 - t153 * (t151 * t168 + t165), -t160 * t148 - t156 * t162, (-t151 * t155 * t156 + t164) * t154 + t156 * t172, t144 * t156 + t145 * t160 + 0; (t147 * t157 + t151 * t171) * t160 + t156 * (-t153 * t155 + t157 * t166), -t156 * t148 + t160 * t162, (t151 * t165 + t168) * t154 - t160 * t172, -t144 * t160 + t145 * t156 + 0; t146 * t157 + t150 * t171, -t146 * t153 + t150 * t167, t151 * t158 + t155 * t170, t150 * t163 - t151 * t161 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end