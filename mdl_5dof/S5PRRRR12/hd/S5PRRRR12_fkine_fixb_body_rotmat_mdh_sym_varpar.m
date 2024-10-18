% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR12 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for the body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRRRR12_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:06:58
	% EndTime: 2024-09-28 18:06:58
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:06:58
	% EndTime: 2024-09-28 18:06:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t99 = cos(pkin(11));
	t98 = sin(pkin(11));
	t1 = [t99, -t98, 0, 0; t98, t99, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:06:58
	% EndTime: 2024-09-28 18:06:58
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t100 = sin(pkin(11));
	t101 = sin(pkin(5));
	t109 = t100 * t101;
	t102 = cos(pkin(11));
	t108 = t102 * t101;
	t103 = cos(pkin(5));
	t104 = sin(qJ(2));
	t107 = t103 * t104;
	t105 = cos(qJ(2));
	t106 = t103 * t105;
	t1 = [-t100 * t107 + t102 * t105, -t100 * t106 - t102 * t104, t109, t102 * pkin(1) + pkin(7) * t109 + 0; t100 * t105 + t102 * t107, -t100 * t104 + t102 * t106, -t108, t100 * pkin(1) - pkin(7) * t108 + 0; t101 * t104, t101 * t105, t103, t103 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:06:58
	% EndTime: 2024-09-28 18:06:58
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (53->22), mult. (37->26), div. (0->0), fcn. (44->12), ass. (0->20)
	t129 = pkin(2) * sin(qJ(2));
	t122 = qJ(2) + qJ(3);
	t128 = pkin(7) + pkin(8);
	t126 = cos(pkin(5));
	t125 = cos(pkin(11));
	t124 = sin(pkin(5));
	t123 = sin(pkin(11));
	t121 = cos(t122);
	t120 = sin(t122);
	t119 = pkin(5) - t122;
	t118 = pkin(5) + t122;
	t117 = cos(qJ(2)) * pkin(2) + pkin(1);
	t116 = cos(t118);
	t115 = sin(t119);
	t114 = cos(t119) / 0.2e1;
	t113 = sin(t118) / 0.2e1;
	t112 = -t124 * t128 + t126 * t129;
	t111 = t114 + t116 / 0.2e1;
	t110 = t113 - t115 / 0.2e1;
	t1 = [-t123 * t110 + t125 * t121, -t123 * t111 - t125 * t120, t123 * t124, -t123 * t112 + t125 * t117 + 0; t125 * t110 + t123 * t121, t125 * t111 - t123 * t120, -t125 * t124, t125 * t112 + t123 * t117 + 0; t114 - t116 / 0.2e1, t113 + t115 / 0.2e1, t126, t124 * t129 + t126 * t128 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:06:58
	% EndTime: 2024-09-28 18:06:58
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (82->26), mult. (51->30), div. (0->0), fcn. (58->15), ass. (0->22)
	t153 = qJ(2) + qJ(3);
	t143 = qJ(4) + t153;
	t151 = cos(qJ(2));
	t152 = pkin(3) * sin(qJ(3)) * t151 + (cos(qJ(3)) * pkin(3) + pkin(2)) * sin(qJ(2));
	t148 = cos(pkin(5));
	t147 = cos(pkin(11));
	t146 = sin(pkin(5));
	t145 = sin(pkin(11));
	t144 = pkin(8) + pkin(9) + pkin(7);
	t141 = cos(t143);
	t140 = sin(t143);
	t139 = pkin(5) - t143;
	t138 = pkin(5) + t143;
	t137 = cos(t138);
	t136 = sin(t139);
	t135 = cos(t139) / 0.2e1;
	t134 = sin(t138) / 0.2e1;
	t133 = pkin(1) + pkin(3) * cos(t153) + t151 * pkin(2);
	t132 = t135 + t137 / 0.2e1;
	t131 = t134 - t136 / 0.2e1;
	t130 = -t146 * t144 + t148 * t152;
	t1 = [-t131 * t145 + t141 * t147, -t132 * t145 - t140 * t147, t145 * t146, -t130 * t145 + t133 * t147 + 0; t131 * t147 + t141 * t145, t132 * t147 - t140 * t145, -t147 * t146, t130 * t147 + t133 * t145 + 0; t135 - t137 / 0.2e1, t134 + t136 / 0.2e1, t148, t148 * t144 + t146 * t152 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:06:58
	% EndTime: 2024-09-28 18:06:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (129->67), mult. (223->115), div. (0->0), fcn. (257->16), ass. (0->46)
	t168 = sin(pkin(6));
	t199 = pkin(10) * t168;
	t167 = sin(pkin(11));
	t172 = cos(pkin(5));
	t198 = t167 * t172;
	t169 = sin(pkin(5));
	t197 = t168 * t169;
	t170 = cos(pkin(11));
	t196 = t168 * t170;
	t195 = t169 * t170;
	t171 = cos(pkin(6));
	t194 = t169 * t171;
	t193 = t170 * t172;
	t173 = sin(qJ(5));
	t192 = t171 * t173;
	t177 = cos(qJ(5));
	t191 = t171 * t177;
	t190 = t172 * t173;
	t189 = t172 * t177;
	t188 = t167 * t199;
	t187 = pkin(10) * t196;
	t186 = t167 * t197;
	t185 = t168 * t195;
	t184 = t171 * t190;
	t183 = t171 * t189;
	t156 = t170 * pkin(4) + t172 * t188;
	t158 = pkin(4) * t198 - t187;
	t174 = sin(qJ(4));
	t178 = cos(qJ(4));
	t182 = t170 * pkin(3) + t156 * t178 - t158 * t174;
	t157 = -t167 * pkin(4) + t172 * t187;
	t159 = pkin(4) * t193 + t188;
	t181 = pkin(3) * t167 - t157 * t178 + t159 * t174;
	t180 = cos(qJ(2));
	t179 = cos(qJ(3));
	t176 = sin(qJ(2));
	t175 = sin(qJ(3));
	t166 = qJ(2) + qJ(3) + qJ(4);
	t165 = cos(t166);
	t164 = sin(t166);
	t162 = t171 * pkin(10) + pkin(7) + pkin(8) + pkin(9);
	t161 = -t174 * pkin(4) + t178 * t199;
	t160 = pkin(4) * t178 + t174 * t199 + pkin(3);
	t155 = pkin(3) * t193 + t157 * t174 + t159 * t178;
	t154 = -pkin(3) * t198 - t156 * t174 - t158 * t178;
	t1 = [(-t167 * t184 + t170 * t177) * t165 + (-t167 * t189 - t170 * t192) * t164 + t173 * t186, (-t167 * t183 - t170 * t173) * t165 + (t167 * t190 - t170 * t191) * t164 + t177 * t186, t164 * t196 + (t165 * t168 * t172 + t194) * t167, (t170 * pkin(2) + t154 * t175 + t182 * t179) * t180 + (-pkin(2) * t198 + t154 * t179 - t182 * t175) * t176 + t169 * t162 * t167 + t170 * pkin(1) + 0; (t167 * t177 + t170 * t184) * t165 + (-t167 * t192 + t170 * t189) * t164 - t173 * t185, (-t167 * t173 + t170 * t183) * t165 + (-t167 * t191 - t170 * t190) * t164 - t177 * t185, -t170 * t194 + (t164 * t167 - t165 * t193) * t168, (t167 * pkin(2) + t155 * t175 + t181 * t179) * t180 + (pkin(2) * t193 + t155 * t179 - t181 * t175) * t176 - t162 * t195 + t167 * pkin(1) + 0; t168 * t190 + (t164 * t177 + t165 * t192) * t169, t168 * t189 + (-t164 * t173 + t165 * t191) * t169, -t165 * t197 + t172 * t171, t162 * t172 + qJ(1) + 0 + ((t160 * t179 + t161 * t175 + pkin(2)) * t176 - (-t175 * t160 + t161 * t179) * t180) * t169; 0, 0, 0, 1;];
	Tc_mdh = t1;
end