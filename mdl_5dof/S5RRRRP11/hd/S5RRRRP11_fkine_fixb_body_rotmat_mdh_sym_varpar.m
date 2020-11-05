% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRP11 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:47
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRRP11_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:47:40
	% EndTime: 2020-11-04 20:47:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:47:40
	% EndTime: 2020-11-04 20:47:40
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t86 = cos(qJ(1));
	t85 = sin(qJ(1));
	t1 = [t86, -t85, 0, 0; t85, t86, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:47:40
	% EndTime: 2020-11-04 20:47:40
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t87 = sin(pkin(5));
	t90 = sin(qJ(1));
	t98 = t90 * t87;
	t89 = sin(qJ(2));
	t97 = t90 * t89;
	t91 = cos(qJ(2));
	t96 = t90 * t91;
	t92 = cos(qJ(1));
	t95 = t92 * t87;
	t94 = t92 * t89;
	t93 = t92 * t91;
	t88 = cos(pkin(5));
	t1 = [-t88 * t97 + t93, -t88 * t96 - t94, t98, t92 * pkin(1) + pkin(7) * t98 + 0; t88 * t94 + t96, t88 * t93 - t97, -t95, t90 * pkin(1) - pkin(7) * t95 + 0; t87 * t89, t87 * t91, t88, t88 * pkin(7) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:47:40
	% EndTime: 2020-11-04 20:47:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t102 = sin(pkin(5));
	t107 = cos(qJ(3));
	t117 = t102 * t107;
	t104 = sin(qJ(3));
	t116 = t104 * t102;
	t105 = sin(qJ(2));
	t115 = t105 * t107;
	t106 = sin(qJ(1));
	t114 = t106 * t105;
	t108 = cos(qJ(2));
	t113 = t106 * t108;
	t109 = cos(qJ(1));
	t112 = t109 * t105;
	t111 = t109 * t108;
	t110 = pkin(2) * t105 - pkin(8) * t108;
	t103 = cos(pkin(5));
	t101 = t108 * pkin(2) + t105 * pkin(8) + pkin(1);
	t100 = t103 * t112 + t113;
	t99 = t102 * pkin(7) - t110 * t103;
	t1 = [(-t103 * t115 + t116) * t106 + t107 * t111, (t103 * t114 - t111) * t104 + t106 * t117, t103 * t113 + t112, t101 * t109 + t99 * t106 + 0; t100 * t107 - t109 * t116, -t100 * t104 - t109 * t117, -t103 * t111 + t114, t101 * t106 - t99 * t109 + 0; t102 * t115 + t103 * t104, t103 * t107 - t105 * t116, -t102 * t108, t103 * pkin(7) + t110 * t102 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:47:40
	% EndTime: 2020-11-04 20:47:40
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (57->34), mult. (112->60), div. (0->0), fcn. (146->10), ass. (0->29)
	t125 = sin(pkin(5));
	t132 = cos(qJ(3));
	t146 = t125 * t132;
	t127 = sin(qJ(4));
	t133 = cos(qJ(2));
	t145 = t127 * t133;
	t128 = sin(qJ(3));
	t144 = t128 * t125;
	t129 = sin(qJ(2));
	t143 = t129 * t132;
	t130 = sin(qJ(1));
	t142 = t130 * t133;
	t131 = cos(qJ(4));
	t141 = t131 * t133;
	t140 = t132 * t133;
	t134 = cos(qJ(1));
	t139 = t134 * t129;
	t138 = t134 * t133;
	t123 = t132 * pkin(3) + t128 * pkin(9) + pkin(2);
	t137 = t133 * pkin(8) - t123 * t129;
	t136 = t128 * pkin(3) - t132 * pkin(9) + pkin(7);
	t126 = cos(pkin(5));
	t121 = t126 * t143 - t144;
	t135 = t121 * t127 + t126 * t141;
	t122 = t127 * t140 - t131 * t129;
	t120 = t125 * t143 + t126 * t128;
	t119 = t129 * pkin(8) + t123 * t133 + pkin(1);
	t118 = t125 * t136 + t137 * t126;
	t1 = [(-t121 * t130 + t132 * t138) * t131 + (t126 * t142 + t139) * t127, -t134 * t122 + t135 * t130, (-t130 * t126 * t129 + t138) * t128 - t130 * t146, t118 * t130 + t119 * t134 + 0; (t121 * t131 - t126 * t145) * t134 + t130 * (t127 * t129 + t131 * t140), -t130 * t122 - t135 * t134, (t126 * t139 + t142) * t128 + t134 * t146, -t118 * t134 + t119 * t130 + 0; t120 * t131 - t125 * t145, -t120 * t127 - t125 * t141, -t126 * t132 + t129 * t144, -t137 * t125 + t136 * t126 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:47:40
	% EndTime: 2020-11-04 20:47:40
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (83->38), mult. (140->64), div. (0->0), fcn. (174->10), ass. (0->31)
	t159 = sin(qJ(4));
	t163 = cos(qJ(4));
	t153 = pkin(4) * t163 + qJ(5) * t159 + pkin(3);
	t160 = sin(qJ(3));
	t164 = cos(qJ(3));
	t180 = -t164 * pkin(9) + t153 * t160 + pkin(7);
	t157 = sin(pkin(5));
	t178 = t157 * t164;
	t165 = cos(qJ(2));
	t177 = t159 * t165;
	t176 = t160 * t157;
	t161 = sin(qJ(2));
	t175 = t161 * t164;
	t162 = sin(qJ(1));
	t174 = t162 * t165;
	t173 = t163 * t165;
	t172 = t164 * t165;
	t166 = cos(qJ(1));
	t171 = t166 * t161;
	t170 = t166 * t165;
	t149 = t160 * pkin(9) + t153 * t164 + pkin(2);
	t154 = -t159 * pkin(4) + qJ(5) * t163 - pkin(8);
	t168 = t149 * t161 + t154 * t165;
	t158 = cos(pkin(5));
	t151 = t158 * t175 - t176;
	t167 = t151 * t159 + t158 * t173;
	t152 = t159 * t172 - t163 * t161;
	t150 = t157 * t175 + t158 * t160;
	t148 = t149 * t165 - t154 * t161 + pkin(1);
	t147 = -t180 * t157 + t168 * t158;
	t1 = [(-t151 * t162 + t164 * t170) * t163 + (t158 * t174 + t171) * t159, (-t162 * t158 * t161 + t170) * t160 - t162 * t178, t166 * t152 - t167 * t162, -t147 * t162 + t148 * t166 + 0; (t151 * t163 - t158 * t177) * t166 + t162 * (t159 * t161 + t163 * t172), (t158 * t171 + t174) * t160 + t166 * t178, t162 * t152 + t167 * t166, t147 * t166 + t148 * t162 + 0; t150 * t163 - t157 * t177, -t158 * t164 + t161 * t176, t150 * t159 + t157 * t173, t168 * t157 + t180 * t158 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end