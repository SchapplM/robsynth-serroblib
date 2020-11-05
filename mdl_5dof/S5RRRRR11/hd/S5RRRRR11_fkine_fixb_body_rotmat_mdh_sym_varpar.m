% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR11 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:51
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRRR11_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR11_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:51:40
	% EndTime: 2020-11-04 20:51:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:51:40
	% EndTime: 2020-11-04 20:51:40
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t80 = cos(qJ(1));
	t79 = sin(qJ(1));
	t1 = [t80, -t79, 0, 0; t79, t80, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:51:40
	% EndTime: 2020-11-04 20:51:40
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t81 = sin(pkin(5));
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
	t82 = cos(pkin(5));
	t1 = [-t82 * t91 + t87, -t82 * t90 - t88, t92, t86 * pkin(1) + pkin(7) * t92 + 0; t82 * t88 + t90, t82 * t87 - t91, -t89, t84 * pkin(1) - pkin(7) * t89 + 0; t81 * t83, t81 * t85, t82, t82 * pkin(7) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:51:40
	% EndTime: 2020-11-04 20:51:40
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (29->23), mult. (60->36), div. (0->0), fcn. (81->8), ass. (0->20)
	t96 = sin(pkin(5));
	t98 = sin(qJ(3));
	t111 = t98 * t96;
	t100 = sin(qJ(1));
	t99 = sin(qJ(2));
	t110 = t100 * t99;
	t103 = cos(qJ(1));
	t109 = t103 * t99;
	t101 = cos(qJ(3));
	t108 = t96 * t101;
	t97 = cos(pkin(5));
	t107 = t97 * t101;
	t102 = cos(qJ(2));
	t106 = t100 * t102;
	t105 = t103 * t102;
	t104 = pkin(2) * t99 - pkin(8) * t102;
	t95 = t102 * pkin(2) + t99 * pkin(8) + pkin(1);
	t94 = t97 * t109 + t106;
	t93 = t96 * pkin(7) - t104 * t97;
	t1 = [(-t99 * t107 + t111) * t100 + t101 * t105, (t97 * t110 - t105) * t98 + t100 * t108, t97 * t106 + t109, t93 * t100 + t95 * t103 + 0; t94 * t101 - t103 * t111, -t103 * t108 - t94 * t98, -t97 * t105 + t110, t95 * t100 - t93 * t103 + 0; t99 * t108 + t97 * t98, -t99 * t111 + t107, -t96 * t102, t97 * pkin(7) + t104 * t96 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:51:40
	% EndTime: 2020-11-04 20:51:40
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (57->34), mult. (112->60), div. (0->0), fcn. (146->10), ass. (0->29)
	t119 = sin(pkin(5));
	t126 = cos(qJ(3));
	t140 = t119 * t126;
	t121 = sin(qJ(4));
	t127 = cos(qJ(2));
	t139 = t121 * t127;
	t122 = sin(qJ(3));
	t138 = t122 * t119;
	t123 = sin(qJ(2));
	t137 = t123 * t126;
	t124 = sin(qJ(1));
	t136 = t124 * t127;
	t125 = cos(qJ(4));
	t135 = t125 * t127;
	t134 = t126 * t127;
	t128 = cos(qJ(1));
	t133 = t128 * t123;
	t132 = t128 * t127;
	t117 = t126 * pkin(3) + t122 * pkin(9) + pkin(2);
	t131 = t127 * pkin(8) - t117 * t123;
	t130 = t122 * pkin(3) - t126 * pkin(9) + pkin(7);
	t120 = cos(pkin(5));
	t115 = t120 * t137 - t138;
	t129 = t115 * t121 + t120 * t135;
	t116 = t121 * t134 - t125 * t123;
	t114 = t119 * t137 + t120 * t122;
	t113 = t123 * pkin(8) + t117 * t127 + pkin(1);
	t112 = t119 * t130 + t131 * t120;
	t1 = [(-t115 * t124 + t126 * t132) * t125 + (t120 * t136 + t133) * t121, -t128 * t116 + t129 * t124, (-t124 * t120 * t123 + t132) * t122 - t124 * t140, t112 * t124 + t113 * t128 + 0; (t115 * t125 - t120 * t139) * t128 + t124 * (t121 * t123 + t125 * t134), -t124 * t116 - t129 * t128, (t120 * t133 + t136) * t122 + t128 * t140, -t112 * t128 + t113 * t124 + 0; t114 * t125 - t119 * t139, -t114 * t121 - t119 * t135, -t120 * t126 + t123 * t138, -t131 * t119 + t130 * t120 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:51:40
	% EndTime: 2020-11-04 20:51:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (90->36), mult. (125->55), div. (0->0), fcn. (159->12), ass. (0->34)
	t154 = sin(pkin(5));
	t159 = cos(qJ(3));
	t174 = t154 * t159;
	t160 = cos(qJ(2));
	t173 = t154 * t160;
	t156 = sin(qJ(3));
	t172 = t156 * t154;
	t157 = sin(qJ(2));
	t171 = t157 * t159;
	t158 = sin(qJ(1));
	t170 = t158 * t157;
	t169 = t158 * t160;
	t161 = cos(qJ(1));
	t168 = t161 * t157;
	t167 = t161 * t160;
	t150 = cos(qJ(4)) * pkin(4) + pkin(3);
	t162 = pkin(10) + pkin(9);
	t143 = t150 * t159 + t162 * t156 + pkin(2);
	t149 = sin(qJ(4)) * pkin(4) + pkin(8);
	t166 = t143 * t157 - t149 * t160;
	t165 = t150 * t156 - t162 * t159 + pkin(7);
	t155 = cos(pkin(5));
	t145 = t155 * t171 - t172;
	t164 = -t145 * t158 + t159 * t167;
	t163 = t145 * t161 + t159 * t169;
	t153 = qJ(4) + qJ(5);
	t152 = cos(t153);
	t151 = sin(t153);
	t147 = t155 * t169 + t168;
	t146 = t155 * t167 - t170;
	t144 = t154 * t171 + t155 * t156;
	t142 = t143 * t160 + t149 * t157 + pkin(1);
	t141 = t165 * t154 - t166 * t155;
	t1 = [t147 * t151 + t164 * t152, t147 * t152 - t164 * t151, (-t155 * t170 + t167) * t156 - t158 * t174, t141 * t158 + t142 * t161 + 0; -t151 * t146 + t163 * t152, -t152 * t146 - t163 * t151, (t155 * t168 + t169) * t156 + t161 * t174, -t141 * t161 + t142 * t158 + 0; t144 * t152 - t151 * t173, -t144 * t151 - t152 * t173, -t155 * t159 + t157 * t172, t166 * t154 + t165 * t155 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end