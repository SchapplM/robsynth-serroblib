% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPR12 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:44
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:28
	% EndTime: 2020-11-04 20:44:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:28
	% EndTime: 2020-11-04 20:44:29
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t81 = cos(qJ(1));
	t80 = sin(qJ(1));
	t1 = [t81, -t80, 0, 0; t80, t81, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:29
	% EndTime: 2020-11-04 20:44:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t82 = sin(pkin(5));
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
	t83 = cos(pkin(5));
	t1 = [-t83 * t92 + t88, -t83 * t91 - t89, t93, t87 * pkin(1) + pkin(7) * t93 + 0; t83 * t89 + t91, t83 * t88 - t92, -t90, t85 * pkin(1) - pkin(7) * t90 + 0; t82 * t84, t82 * t86, t83, t83 * pkin(7) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:29
	% EndTime: 2020-11-04 20:44:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t97 = sin(pkin(5));
	t99 = sin(qJ(3));
	t112 = t99 * t97;
	t102 = cos(qJ(3));
	t111 = t97 * t102;
	t100 = sin(qJ(2));
	t110 = t100 * t102;
	t101 = sin(qJ(1));
	t109 = t101 * t100;
	t103 = cos(qJ(2));
	t108 = t101 * t103;
	t104 = cos(qJ(1));
	t107 = t104 * t100;
	t106 = t104 * t103;
	t105 = pkin(2) * t100 - pkin(8) * t103;
	t98 = cos(pkin(5));
	t96 = t103 * pkin(2) + t100 * pkin(8) + pkin(1);
	t95 = t98 * t107 + t108;
	t94 = t97 * pkin(7) - t105 * t98;
	t1 = [(-t98 * t110 + t112) * t101 + t102 * t106, (t98 * t109 - t106) * t99 + t101 * t111, t98 * t108 + t107, t94 * t101 + t96 * t104 + 0; t95 * t102 - t104 * t112, -t104 * t111 - t95 * t99, -t98 * t106 + t109, t96 * t101 - t94 * t104 + 0; t97 * t110 + t98 * t99, -t100 * t112 + t98 * t102, -t97 * t103, t98 * pkin(7) + t105 * t97 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:29
	% EndTime: 2020-11-04 20:44:29
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (57->32), mult. (112->56), div. (0->0), fcn. (146->10), ass. (0->29)
	t123 = sin(pkin(5));
	t129 = cos(qJ(3));
	t140 = t123 * t129;
	t130 = cos(qJ(2));
	t139 = t123 * t130;
	t125 = cos(pkin(5));
	t127 = sin(qJ(2));
	t138 = t125 * t127;
	t137 = t125 * t130;
	t126 = sin(qJ(3));
	t136 = t126 * t123;
	t135 = t127 * t129;
	t134 = t129 * t130;
	t121 = pkin(3) * t129 + qJ(4) * t126 + pkin(2);
	t133 = t130 * pkin(8) - t121 * t127;
	t132 = t126 * pkin(3) - qJ(4) * t129 + pkin(7);
	t131 = cos(qJ(1));
	t128 = sin(qJ(1));
	t124 = cos(pkin(10));
	t122 = sin(pkin(10));
	t120 = t122 * t127 + t124 * t134;
	t119 = t125 * t135 - t136;
	t118 = t123 * t135 + t125 * t126;
	t117 = t122 * t134 - t127 * t124;
	t116 = t127 * pkin(8) + t121 * t130 + pkin(1);
	t115 = t122 * t119 + t124 * t137;
	t114 = -t124 * t119 + t122 * t137;
	t113 = t123 * t132 + t133 * t125;
	t1 = [t114 * t128 + t131 * t120, t115 * t128 - t131 * t117, (-t128 * t138 + t131 * t130) * t126 - t128 * t140, t113 * t128 + t116 * t131 + 0; -t114 * t131 + t128 * t120, -t115 * t131 - t128 * t117, (t128 * t130 + t131 * t138) * t126 + t131 * t140, -t113 * t131 + t116 * t128 + 0; t118 * t124 - t122 * t139, -t118 * t122 - t124 * t139, -t125 * t129 + t127 * t136, -t133 * t123 + t132 * t125 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:29
	% EndTime: 2020-11-04 20:44:29
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (90->36), mult. (125->55), div. (0->0), fcn. (159->12), ass. (0->34)
	t154 = sin(pkin(5));
	t160 = cos(qJ(3));
	t174 = t154 * t160;
	t161 = cos(qJ(2));
	t173 = t154 * t161;
	t157 = sin(qJ(3));
	t172 = t157 * t154;
	t158 = sin(qJ(2));
	t171 = t158 * t160;
	t159 = sin(qJ(1));
	t170 = t159 * t158;
	t169 = t159 * t161;
	t162 = cos(qJ(1));
	t168 = t162 * t158;
	t167 = t162 * t161;
	t150 = cos(pkin(10)) * pkin(4) + pkin(3);
	t156 = qJ(4) + pkin(9);
	t143 = t150 * t160 + t156 * t157 + pkin(2);
	t149 = sin(pkin(10)) * pkin(4) + pkin(8);
	t166 = t143 * t158 - t149 * t161;
	t165 = t150 * t157 - t156 * t160 + pkin(7);
	t155 = cos(pkin(5));
	t145 = t155 * t171 - t172;
	t164 = -t145 * t159 + t160 * t167;
	t163 = t145 * t162 + t160 * t169;
	t153 = pkin(10) + qJ(5);
	t152 = cos(t153);
	t151 = sin(t153);
	t147 = t155 * t169 + t168;
	t146 = t155 * t167 - t170;
	t144 = t154 * t171 + t155 * t157;
	t142 = t143 * t161 + t149 * t158 + pkin(1);
	t141 = t165 * t154 - t166 * t155;
	t1 = [t147 * t151 + t164 * t152, t147 * t152 - t164 * t151, (-t155 * t170 + t167) * t157 - t159 * t174, t141 * t159 + t142 * t162 + 0; -t151 * t146 + t163 * t152, -t152 * t146 - t163 * t151, (t155 * t168 + t169) * t157 + t162 * t174, -t141 * t162 + t142 * t159 + 0; t144 * t152 - t151 * t173, -t144 * t151 - t152 * t173, -t155 * t160 + t158 * t172, t166 * t154 + t165 * t155 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end