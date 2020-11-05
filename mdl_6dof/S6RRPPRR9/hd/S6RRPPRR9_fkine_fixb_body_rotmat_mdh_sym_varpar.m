% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR9 (for one body)
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
% Datum: 2020-11-04 22:05
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:09
	% EndTime: 2020-11-04 22:05:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:09
	% EndTime: 2020-11-04 22:05:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t80 = cos(qJ(1));
	t79 = sin(qJ(1));
	t1 = [t80, -t79, 0, 0; t79, t80, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:09
	% EndTime: 2020-11-04 22:05:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t81 = sin(pkin(6));
	t84 = sin(qJ(1));
	t92 = t84 * t81;
	t83 = sin(qJ(2));
	t91 = t84 * t83;
	t85 = cos(qJ(2));
	t90 = t84 * t85;
	t86 = cos(qJ(1));
	t89 = t86 * t81;
	t88 = t86 * t83;
	t87 = t86 * t85;
	t82 = cos(pkin(6));
	t1 = [-t82 * t91 + t87, -t82 * t90 - t88, t92, t86 * pkin(1) + pkin(8) * t92 + 0; t82 * t88 + t90, t82 * t87 - t91, -t89, t84 * pkin(1) - pkin(8) * t89 + 0; t81 * t83, t81 * t85, t82, t82 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:09
	% EndTime: 2020-11-04 22:05:10
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (23->18), mult. (39->24), div. (0->0), fcn. (52->6), ass. (0->14)
	t97 = sin(qJ(2));
	t98 = sin(qJ(1));
	t105 = t98 * t97;
	t99 = cos(qJ(2));
	t104 = t98 * t99;
	t100 = cos(qJ(1));
	t103 = t100 * t97;
	t102 = t100 * t99;
	t101 = pkin(2) * t97 - qJ(3) * t99;
	t96 = cos(pkin(6));
	t95 = sin(pkin(6));
	t94 = pkin(2) * t99 + qJ(3) * t97 + pkin(1);
	t93 = t95 * pkin(8) - t101 * t96;
	t1 = [t98 * t95, t105 * t96 - t102, t104 * t96 + t103, t100 * t94 + t93 * t98 + 0; -t100 * t95, -t103 * t96 - t104, -t102 * t96 + t105, -t100 * t93 + t94 * t98 + 0; t96, -t95 * t97, -t95 * t99, t96 * pkin(8) + t101 * t95 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:10
	% EndTime: 2020-11-04 22:05:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->19), mult. (39->24), div. (0->0), fcn. (52->6), ass. (0->16)
	t111 = sin(qJ(2));
	t112 = sin(qJ(1));
	t120 = t112 * t111;
	t113 = cos(qJ(2));
	t119 = t112 * t113;
	t114 = cos(qJ(1));
	t118 = t114 * t111;
	t117 = t114 * t113;
	t110 = pkin(2) + qJ(4);
	t116 = qJ(3) * t113 - t110 * t111;
	t115 = pkin(3) + pkin(8);
	t109 = cos(pkin(6));
	t108 = sin(pkin(6));
	t107 = t111 * qJ(3) + t110 * t113 + pkin(1);
	t106 = t108 * t115 + t116 * t109;
	t1 = [t112 * t108, t109 * t119 + t118, -t109 * t120 + t117, t106 * t112 + t107 * t114 + 0; -t114 * t108, -t109 * t117 + t120, t109 * t118 + t119, -t106 * t114 + t107 * t112 + 0; t109, -t108 * t113, t108 * t111, -t116 * t108 + t115 * t109 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:10
	% EndTime: 2020-11-04 22:05:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (44->25), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->23)
	t125 = sin(pkin(6));
	t132 = cos(qJ(5));
	t142 = t125 * t132;
	t129 = sin(qJ(5));
	t141 = t129 * t125;
	t130 = sin(qJ(2));
	t140 = t130 * t132;
	t131 = sin(qJ(1));
	t139 = t131 * t130;
	t133 = cos(qJ(2));
	t138 = t131 * t133;
	t134 = cos(qJ(1));
	t137 = t134 * t130;
	t136 = t134 * t133;
	t127 = qJ(3) - pkin(9);
	t128 = pkin(2) + qJ(4);
	t135 = t127 * t133 - t128 * t130;
	t126 = cos(pkin(6));
	t124 = pkin(3) + pkin(4) + pkin(8);
	t123 = t127 * t130 + t128 * t133 + pkin(1);
	t122 = t126 * t137 + t138;
	t121 = t125 * t124 + t135 * t126;
	t1 = [t131 * t142 + (-t126 * t139 + t136) * t129, (-t126 * t140 - t141) * t131 + t132 * t136, -t126 * t138 - t137, t121 * t131 + t123 * t134 + 0; t122 * t129 - t134 * t142, t122 * t132 + t134 * t141, t126 * t136 - t139, -t121 * t134 + t123 * t131 + 0; t126 * t132 + t130 * t141, t125 * t140 - t126 * t129, t125 * t133, t124 * t126 - t135 * t125 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:10
	% EndTime: 2020-11-04 22:05:10
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (73->38), mult. (112->58), div. (0->0), fcn. (146->10), ass. (0->30)
	t151 = cos(pkin(6));
	t158 = cos(qJ(5));
	t172 = t151 * t158;
	t153 = sin(qJ(6));
	t159 = cos(qJ(2));
	t171 = t153 * t159;
	t150 = sin(pkin(6));
	t154 = sin(qJ(5));
	t170 = t154 * t150;
	t155 = sin(qJ(2));
	t169 = t154 * t155;
	t156 = sin(qJ(1));
	t168 = t156 * t159;
	t157 = cos(qJ(6));
	t167 = t157 * t159;
	t166 = t158 * t150;
	t160 = cos(qJ(1));
	t165 = t159 * t160;
	t164 = t160 * t155;
	t148 = t154 * pkin(5) - pkin(10) * t158 + pkin(2) + qJ(4);
	t152 = qJ(3) - pkin(9);
	t163 = t148 * t155 - t152 * t159;
	t145 = t151 * t169 - t166;
	t162 = -t145 * t153 + t151 * t167;
	t161 = t158 * pkin(5) + t154 * pkin(10) + pkin(3) + pkin(4) + pkin(8);
	t147 = t154 * t171 + t157 * t155;
	t146 = t150 * t169 + t172;
	t144 = t148 * t159 + t152 * t155 + pkin(1);
	t143 = t150 * t161 - t163 * t151;
	t1 = [(-t145 * t156 + t154 * t165) * t157 - t153 * (t151 * t168 + t164), -t160 * t147 - t162 * t156, (t155 * t172 + t170) * t156 - t158 * t165, t143 * t156 + t144 * t160 + 0; (t145 * t157 + t151 * t171) * t160 + t156 * (-t153 * t155 + t154 * t167), -t156 * t147 + t162 * t160, (-t151 * t164 - t168) * t158 - t160 * t170, -t143 * t160 + t144 * t156 + 0; t146 * t157 + t150 * t171, -t146 * t153 + t150 * t167, t151 * t154 - t155 * t166, t163 * t150 + t161 * t151 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end