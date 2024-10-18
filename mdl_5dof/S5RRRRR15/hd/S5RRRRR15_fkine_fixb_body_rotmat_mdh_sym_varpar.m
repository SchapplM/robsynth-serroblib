% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR15 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for the body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRRR15_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR15_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:23:44
	% EndTime: 2024-09-27 22:23:44
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:23:44
	% EndTime: 2024-09-27 22:23:44
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t96 = cos(qJ(1));
	t95 = sin(qJ(1));
	t1 = [t96, -t95, 0, 0; t95, t96, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:23:44
	% EndTime: 2024-09-27 22:23:44
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t100 = sin(qJ(1));
	t97 = sin(pkin(5));
	t108 = t100 * t97;
	t99 = sin(qJ(2));
	t107 = t100 * t99;
	t102 = cos(qJ(1));
	t106 = t102 * t97;
	t105 = t102 * t99;
	t101 = cos(qJ(2));
	t104 = t100 * t101;
	t103 = t102 * t101;
	t98 = cos(pkin(5));
	t1 = [-t98 * t107 + t103, -t98 * t104 - t105, t108, t102 * pkin(1) + pkin(8) * t108 + 0; t98 * t105 + t104, t98 * t103 - t107, -t106, t100 * pkin(1) - pkin(8) * t106 + 0; t97 * t99, t97 * t101, t98, t98 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:23:44
	% EndTime: 2024-09-27 22:23:44
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (53->22), mult. (33->26), div. (0->0), fcn. (44->12), ass. (0->22)
	t125 = sin(qJ(1));
	t130 = -t125 / 0.2e1;
	t126 = cos(qJ(1));
	t129 = t126 / 0.2e1;
	t128 = pkin(2) * sin(qJ(2));
	t121 = qJ(2) + qJ(3);
	t127 = pkin(8) + pkin(9);
	t123 = cos(pkin(5));
	t122 = sin(pkin(5));
	t120 = cos(t121);
	t119 = sin(t121);
	t118 = pkin(5) - t121;
	t117 = pkin(5) + t121;
	t116 = cos(qJ(2)) * pkin(2) + pkin(1);
	t115 = cos(t118);
	t114 = cos(t117);
	t113 = sin(t118);
	t112 = sin(t117);
	t111 = t115 + t114;
	t110 = t112 - t113;
	t109 = -t122 * t127 + t123 * t128;
	t1 = [t110 * t130 + t126 * t120, t111 * t130 - t126 * t119, t125 * t122, -t109 * t125 + t126 * t116 + 0; t110 * t129 + t125 * t120, t111 * t129 - t125 * t119, -t126 * t122, t126 * t109 + t125 * t116 + 0; t115 / 0.2e1 - t114 / 0.2e1, t112 / 0.2e1 + t113 / 0.2e1, t123, t122 * t128 + t123 * t127 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:23:44
	% EndTime: 2024-09-27 22:23:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (82->26), mult. (51->30), div. (0->0), fcn. (58->15), ass. (0->22)
	t154 = qJ(2) + qJ(3);
	t144 = qJ(4) + t154;
	t151 = cos(qJ(2));
	t153 = pkin(3) * sin(qJ(3)) * t151 + (cos(qJ(3)) * pkin(3) + pkin(2)) * sin(qJ(2));
	t152 = cos(qJ(1));
	t150 = sin(qJ(1));
	t147 = cos(pkin(5));
	t146 = sin(pkin(5));
	t145 = pkin(8) + pkin(9) + pkin(10);
	t142 = cos(t144);
	t141 = sin(t144);
	t140 = pkin(5) - t144;
	t139 = pkin(5) + t144;
	t138 = cos(t139);
	t137 = sin(t140);
	t136 = cos(t140) / 0.2e1;
	t135 = sin(t139) / 0.2e1;
	t134 = pkin(1) + pkin(3) * cos(t154) + t151 * pkin(2);
	t133 = t136 + t138 / 0.2e1;
	t132 = t135 - t137 / 0.2e1;
	t131 = -t146 * t145 + t153 * t147;
	t1 = [-t150 * t132 + t152 * t142, -t150 * t133 - t152 * t141, t150 * t146, -t150 * t131 + t152 * t134 + 0; t152 * t132 + t150 * t142, t152 * t133 - t150 * t141, -t152 * t146, t152 * t131 + t150 * t134 + 0; t136 - t138 / 0.2e1, t135 + t137 / 0.2e1, t147, t147 * t145 + t153 * t146 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:23:44
	% EndTime: 2024-09-27 22:23:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (129->46), mult. (189->76), div. (0->0), fcn. (223->16), ass. (0->42)
	t165 = sin(pkin(6));
	t195 = pkin(11) * t165;
	t164 = qJ(2) + qJ(3) + qJ(4);
	t163 = cos(t164);
	t167 = cos(pkin(6));
	t194 = t163 * t167;
	t168 = cos(pkin(5));
	t193 = t163 * t168;
	t192 = t165 * t168;
	t166 = sin(pkin(5));
	t173 = sin(qJ(1));
	t191 = t166 * t173;
	t178 = cos(qJ(1));
	t190 = t166 * t178;
	t169 = sin(qJ(5));
	t189 = t173 * t169;
	t174 = cos(qJ(5));
	t188 = t174 * t173;
	t187 = t178 * t169;
	t186 = t178 * t174;
	t185 = t168 * t187;
	t184 = t165 * t191;
	t183 = t168 * t189;
	t182 = t168 * t188;
	t181 = t165 * t190;
	t180 = t168 * t186;
	t170 = sin(qJ(4));
	t175 = cos(qJ(4));
	t159 = pkin(4) * t175 + t170 * t195 + pkin(3);
	t160 = -pkin(4) * t170 + t175 * t195;
	t171 = sin(qJ(3));
	t176 = cos(qJ(3));
	t157 = t159 * t176 + t160 * t171 + pkin(2);
	t158 = -t159 * t171 + t160 * t176;
	t172 = sin(qJ(2));
	t177 = cos(qJ(2));
	t179 = t157 * t172 - t158 * t177;
	t162 = sin(t164);
	t161 = pkin(11) * t167 + pkin(8) + pkin(9) + pkin(10);
	t156 = t157 * t177 + t158 * t172 + pkin(1);
	t155 = t166 * t161 - t168 * t179;
	t1 = [(-t167 * t183 + t186) * t163 + (-t167 * t187 - t182) * t162 + t169 * t184, (-t167 * t182 - t187) * t163 + (-t167 * t186 + t183) * t162 + t174 * t184, t167 * t191 + (t162 * t178 + t173 * t193) * t165, t155 * t173 + t156 * t178 + 0; (t167 * t185 + t188) * t163 + (-t167 * t189 + t180) * t162 - t169 * t181, (t167 * t180 - t189) * t163 + (-t167 * t188 - t185) * t162 - t174 * t181, -t167 * t190 + (t162 * t173 - t178 * t193) * t165, -t155 * t178 + t156 * t173 + 0; t169 * t192 + (t162 * t174 + t169 * t194) * t166, t174 * t192 + (-t162 * t169 + t174 * t194) * t166, -t163 * t165 * t166 + t167 * t168, t161 * t168 + t166 * t179 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end