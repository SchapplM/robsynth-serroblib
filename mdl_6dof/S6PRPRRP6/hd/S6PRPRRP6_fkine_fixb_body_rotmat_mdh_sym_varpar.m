% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRP6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:03
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:03:17
	% EndTime: 2020-11-04 21:03:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:03:17
	% EndTime: 2020-11-04 21:03:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t104 = cos(pkin(10));
	t103 = sin(pkin(10));
	t1 = [t104, -t103, 0, 0; t103, t104, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:03:17
	% EndTime: 2020-11-04 21:03:18
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t105 = sin(pkin(10));
	t106 = sin(pkin(6));
	t114 = t105 * t106;
	t107 = cos(pkin(10));
	t113 = t107 * t106;
	t108 = cos(pkin(6));
	t109 = sin(qJ(2));
	t112 = t108 * t109;
	t110 = cos(qJ(2));
	t111 = t108 * t110;
	t1 = [-t105 * t112 + t107 * t110, -t105 * t111 - t107 * t109, t114, t107 * pkin(1) + pkin(7) * t114 + 0; t105 * t110 + t107 * t112, -t105 * t109 + t107 * t111, -t113, t105 * pkin(1) - pkin(7) * t113 + 0; t106 * t109, t106 * t110, t108, t108 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:03:18
	% EndTime: 2020-11-04 21:03:18
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (23->23), mult. (43->36), div. (0->0), fcn. (56->6), ass. (0->13)
	t115 = sin(pkin(10));
	t116 = sin(pkin(6));
	t126 = t115 * t116;
	t118 = cos(pkin(6));
	t125 = t115 * t118;
	t117 = cos(pkin(10));
	t124 = t117 * t116;
	t123 = t117 * t118;
	t119 = sin(qJ(2));
	t122 = t118 * t119;
	t120 = cos(qJ(2));
	t121 = t118 * t120;
	t1 = [t126, t115 * t122 - t117 * t120, t115 * t121 + t117 * t119, (t117 * pkin(2) + qJ(3) * t125) * t120 + (-pkin(2) * t125 + t117 * qJ(3)) * t119 + pkin(7) * t126 + t117 * pkin(1) + 0; -t124, -t115 * t120 - t117 * t122, t115 * t119 - t117 * t121, (t115 * pkin(2) - qJ(3) * t123) * t120 + (pkin(2) * t123 + t115 * qJ(3)) * t119 - pkin(7) * t124 + t115 * pkin(1) + 0; t118, -t116 * t119, -t116 * t120, t118 * pkin(7) + qJ(1) + 0 + (pkin(2) * t119 - qJ(3) * t120) * t116; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:03:18
	% EndTime: 2020-11-04 21:03:18
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (36->28), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->22)
	t129 = sin(pkin(10));
	t147 = t129 * qJ(3);
	t130 = sin(pkin(6));
	t133 = sin(qJ(4));
	t146 = t130 * t133;
	t135 = cos(qJ(4));
	t145 = t130 * t135;
	t136 = cos(qJ(2));
	t144 = t130 * t136;
	t137 = pkin(3) + pkin(7);
	t143 = t130 * t137;
	t131 = cos(pkin(10));
	t142 = t131 * qJ(3);
	t132 = cos(pkin(6));
	t134 = sin(qJ(2));
	t141 = t132 * t134;
	t140 = t132 * t136;
	t138 = pkin(2) + pkin(8);
	t139 = t132 * t138;
	t128 = t129 * t140 + t131 * t134;
	t127 = t129 * t134 - t131 * t140;
	t1 = [t128 * t133 + t129 * t145, t128 * t135 - t129 * t146, -t129 * t141 + t131 * t136, (t131 * t138 + t132 * t147) * t136 + (-t129 * t139 + t142) * t134 + t129 * t143 + t131 * pkin(1) + 0; t127 * t133 - t131 * t145, t127 * t135 + t131 * t146, t129 * t136 + t131 * t141, (t129 * t138 - t132 * t142) * t136 + (t131 * t139 + t147) * t134 - t131 * t143 + t129 * pkin(1) + 0; t132 * t135 - t133 * t144, -t132 * t133 - t135 * t144, t130 * t134, t137 * t132 + qJ(1) + 0 + (-qJ(3) * t136 + t134 * t138) * t130; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:03:18
	% EndTime: 2020-11-04 21:03:18
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (65->34), mult. (140->56), div. (0->0), fcn. (174->10), ass. (0->30)
	t154 = sin(pkin(6));
	t156 = cos(pkin(6));
	t159 = sin(qJ(2));
	t162 = cos(qJ(2));
	t164 = pkin(2) + pkin(8);
	t158 = sin(qJ(4));
	t161 = cos(qJ(4));
	t169 = t158 * pkin(4) - pkin(9) * t161 + qJ(3);
	t166 = -t164 * t159 + t169 * t162;
	t168 = t161 * pkin(4) + t158 * pkin(9) + pkin(3) + pkin(7);
	t180 = t168 * t154 + t166 * t156;
	t177 = t154 * t158;
	t176 = t154 * t159;
	t175 = t156 * t159;
	t174 = t156 * t162;
	t173 = t158 * t159;
	t172 = t158 * t162;
	t171 = t161 * t154;
	t167 = t156 * t172 + t171;
	t165 = t169 * t159 + t164 * t162 + pkin(1);
	t160 = cos(qJ(5));
	t157 = sin(qJ(5));
	t155 = cos(pkin(10));
	t153 = sin(pkin(10));
	t152 = -t154 * t172 + t156 * t161;
	t151 = t153 * t175 - t155 * t162;
	t150 = t153 * t162 + t155 * t175;
	t149 = -t153 * t173 + t167 * t155;
	t148 = t167 * t153 + t155 * t173;
	t1 = [t148 * t160 - t157 * t151, -t148 * t157 - t160 * t151, t153 * t177 - (t153 * t174 + t155 * t159) * t161, t180 * t153 + t165 * t155 + 0; -t149 * t160 + t150 * t157, t149 * t157 + t150 * t160, -t155 * t177 + (-t153 * t159 + t155 * t174) * t161, t165 * t153 - t180 * t155 + 0; t152 * t160 + t157 * t176, -t152 * t157 + t160 * t176, t156 * t158 + t162 * t171, -t166 * t154 + t168 * t156 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:03:18
	% EndTime: 2020-11-04 21:03:18
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (91->38), mult. (178->60), div. (0->0), fcn. (212->10), ass. (0->31)
	t188 = sin(pkin(6));
	t190 = cos(pkin(6));
	t191 = sin(qJ(5));
	t194 = cos(qJ(5));
	t186 = pkin(5) * t194 + qJ(6) * t191 + pkin(4);
	t192 = sin(qJ(4));
	t195 = cos(qJ(4));
	t202 = t192 * pkin(9) + t186 * t195 + pkin(3) + pkin(7);
	t193 = sin(qJ(2));
	t196 = cos(qJ(2));
	t203 = pkin(9) * t195 - t186 * t192 - qJ(3);
	t213 = -t191 * pkin(5) + qJ(6) * t194 - pkin(2) - pkin(8);
	t216 = -t213 * t193 + t203 * t196;
	t217 = -t202 * t188 + t216 * t190;
	t212 = t188 * t192;
	t211 = t188 * t193;
	t210 = t188 * t195;
	t209 = t190 * t193;
	t208 = t190 * t196;
	t207 = t192 * t193;
	t206 = t192 * t196;
	t201 = t190 * t206 + t210;
	t199 = -t203 * t193 - t196 * t213 + pkin(1);
	t189 = cos(pkin(10));
	t187 = sin(pkin(10));
	t185 = -t188 * t206 + t190 * t195;
	t184 = t187 * t209 - t189 * t196;
	t183 = t187 * t196 + t189 * t209;
	t182 = t187 * t207 - t201 * t189;
	t181 = t201 * t187 + t189 * t207;
	t1 = [t181 * t194 - t191 * t184, t187 * t212 - (t187 * t208 + t189 * t193) * t195, t181 * t191 + t194 * t184, -t217 * t187 + t199 * t189 + 0; t182 * t194 + t183 * t191, -t189 * t212 + (-t187 * t193 + t189 * t208) * t195, t182 * t191 - t183 * t194, t199 * t187 + t217 * t189 + 0; t185 * t194 + t191 * t211, t190 * t192 + t196 * t210, t185 * t191 - t194 * t211, t216 * t188 + t202 * t190 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end