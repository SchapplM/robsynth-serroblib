% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:11
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:51
	% EndTime: 2020-11-04 21:11:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:51
	% EndTime: 2020-11-04 21:11:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t94 = cos(pkin(11));
	t93 = sin(pkin(11));
	t1 = [t94, -t93, 0, 0; t93, t94, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:51
	% EndTime: 2020-11-04 21:11:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->18), div. (0->0), fcn. (36->6), ass. (0->12)
	t95 = sin(pkin(11));
	t96 = sin(pkin(6));
	t105 = t95 * t96;
	t97 = cos(pkin(11));
	t104 = t97 * t96;
	t98 = cos(pkin(6));
	t99 = sin(qJ(2));
	t103 = t98 * t99;
	t100 = cos(qJ(2));
	t102 = t95 * t100;
	t101 = t97 * t100;
	t1 = [-t95 * t103 + t101, -t98 * t102 - t97 * t99, t105, t97 * pkin(1) + pkin(7) * t105 + 0; t97 * t103 + t102, t98 * t101 - t95 * t99, -t104, t95 * pkin(1) - pkin(7) * t104 + 0; t96 * t99, t96 * t100, t98, t98 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:51
	% EndTime: 2020-11-04 21:11:51
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t109 = sin(pkin(6));
	t122 = t109 * pkin(7);
	t108 = sin(pkin(11));
	t111 = cos(pkin(6));
	t121 = t108 * t111;
	t112 = sin(qJ(3));
	t120 = t109 * t112;
	t114 = cos(qJ(3));
	t119 = t109 * t114;
	t110 = cos(pkin(11));
	t118 = t110 * t111;
	t113 = sin(qJ(2));
	t117 = t111 * t113;
	t115 = cos(qJ(2));
	t116 = t111 * t115;
	t107 = t108 * t115 + t110 * t117;
	t106 = t108 * t117 - t110 * t115;
	t1 = [-t106 * t114 + t108 * t120, t106 * t112 + t108 * t119, t108 * t116 + t110 * t113, (t110 * pkin(2) + pkin(8) * t121) * t115 + (-pkin(2) * t121 + t110 * pkin(8)) * t113 + t108 * t122 + t110 * pkin(1) + 0; t107 * t114 - t110 * t120, -t107 * t112 - t110 * t119, t108 * t113 - t110 * t116, (t108 * pkin(2) - pkin(8) * t118) * t115 + (pkin(2) * t118 + t108 * pkin(8)) * t113 - t110 * t122 + t108 * pkin(1) + 0; t111 * t112 + t113 * t119, t111 * t114 - t113 * t120, -t109 * t115, t111 * pkin(7) + qJ(1) + 0 + (pkin(2) * t113 - pkin(8) * t115) * t109; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:51
	% EndTime: 2020-11-04 21:11:51
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (57->39), mult. (134->66), div. (0->0), fcn. (178->10), ass. (0->29)
	t135 = sin(pkin(6));
	t150 = t135 * pkin(7);
	t134 = sin(pkin(11));
	t138 = cos(pkin(6));
	t149 = t134 * t138;
	t139 = sin(qJ(3));
	t148 = t135 * t139;
	t141 = cos(qJ(3));
	t147 = t135 * t141;
	t142 = cos(qJ(2));
	t146 = t135 * t142;
	t137 = cos(pkin(11));
	t145 = t137 * t138;
	t140 = sin(qJ(2));
	t144 = t138 * t140;
	t143 = t138 * t142;
	t136 = cos(pkin(12));
	t133 = sin(pkin(12));
	t132 = t138 * t139 + t140 * t147;
	t131 = -t138 * t141 + t140 * t148;
	t130 = t134 * t143 + t137 * t140;
	t129 = t134 * t142 + t137 * t144;
	t128 = t134 * t140 - t137 * t143;
	t127 = t134 * t144 - t137 * t142;
	t126 = -t127 * t141 + t134 * t148;
	t125 = t129 * t141 - t137 * t148;
	t124 = t129 * t139 + t137 * t147;
	t123 = t127 * t139 + t134 * t147;
	t1 = [t126 * t136 + t130 * t133, -t126 * t133 + t130 * t136, -t123, t126 * pkin(3) - t123 * qJ(4) + (t137 * pkin(2) + pkin(8) * t149) * t142 + (-pkin(2) * t149 + t137 * pkin(8)) * t140 + t134 * t150 + t137 * pkin(1) + 0; t125 * t136 + t128 * t133, -t125 * t133 + t128 * t136, t124, t125 * pkin(3) + t124 * qJ(4) + (t134 * pkin(2) - pkin(8) * t145) * t142 + (pkin(2) * t145 + t134 * pkin(8)) * t140 - t137 * t150 + t134 * pkin(1) + 0; t132 * t136 - t133 * t146, -t132 * t133 - t136 * t146, t131, t132 * pkin(3) + t138 * pkin(7) + t131 * qJ(4) + qJ(1) + 0 + (pkin(2) * t140 - pkin(8) * t142) * t135; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:51
	% EndTime: 2020-11-04 21:11:51
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (90->36), mult. (149->58), div. (0->0), fcn. (183->12), ass. (0->33)
	t163 = sin(pkin(6));
	t165 = cos(pkin(6));
	t157 = sin(pkin(12)) * pkin(4) + pkin(8);
	t168 = sin(qJ(2));
	t170 = cos(qJ(2));
	t158 = cos(pkin(12)) * pkin(4) + pkin(3);
	t166 = qJ(4) + pkin(9);
	t167 = sin(qJ(3));
	t169 = cos(qJ(3));
	t174 = t158 * t169 + t166 * t167 + pkin(2);
	t172 = -t157 * t170 + t174 * t168;
	t173 = t158 * t167 - t166 * t169 + pkin(7);
	t185 = t173 * t163 - t172 * t165;
	t181 = t163 * t167;
	t180 = t163 * t169;
	t179 = t163 * t170;
	t178 = t165 * t168;
	t177 = t165 * t170;
	t176 = t168 * t169;
	t175 = t169 * t170;
	t171 = t157 * t168 + t174 * t170 + pkin(1);
	t164 = cos(pkin(11));
	t162 = sin(pkin(11));
	t161 = pkin(12) + qJ(5);
	t160 = cos(t161);
	t159 = sin(t161);
	t156 = t163 * t176 + t165 * t167;
	t155 = -t165 * t176 + t181;
	t154 = t162 * t177 + t164 * t168;
	t153 = t162 * t168 - t164 * t177;
	t152 = t162 * t155 + t164 * t175;
	t151 = -t164 * t155 + t162 * t175;
	t1 = [t152 * t160 + t154 * t159, -t152 * t159 + t154 * t160, -(t162 * t178 - t164 * t170) * t167 - t162 * t180, t185 * t162 + t171 * t164 + 0; t151 * t160 + t153 * t159, -t151 * t159 + t153 * t160, (t162 * t170 + t164 * t178) * t167 + t164 * t180, t171 * t162 - t185 * t164 + 0; t156 * t160 - t159 * t179, -t156 * t159 - t160 * t179, -t165 * t169 + t168 * t181, t172 * t163 + t173 * t165 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:51
	% EndTime: 2020-11-04 21:11:51
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (107->48), mult. (156->72), div. (0->0), fcn. (202->14), ass. (0->34)
	t204 = sin(pkin(6));
	t218 = t204 * pkin(7);
	t203 = sin(pkin(11));
	t206 = cos(pkin(6));
	t217 = t203 * t206;
	t207 = sin(qJ(3));
	t216 = t204 * t207;
	t209 = cos(qJ(3));
	t215 = t204 * t209;
	t210 = cos(qJ(2));
	t214 = t204 * t210;
	t205 = cos(pkin(11));
	t213 = t205 * t206;
	t208 = sin(qJ(2));
	t212 = t206 * t208;
	t211 = t206 * t210;
	t202 = pkin(12) + qJ(5);
	t201 = -pkin(10) - pkin(9) - qJ(4);
	t200 = qJ(6) + t202;
	t199 = cos(t200);
	t198 = sin(t200);
	t197 = pkin(5) * sin(t202) + sin(pkin(12)) * pkin(4);
	t196 = pkin(5) * cos(t202) + cos(pkin(12)) * pkin(4) + pkin(3);
	t195 = t206 * t207 + t208 * t215;
	t194 = -t206 * t209 + t208 * t216;
	t193 = t203 * t211 + t205 * t208;
	t192 = t203 * t210 + t205 * t212;
	t191 = t203 * t208 - t205 * t211;
	t190 = t203 * t212 - t205 * t210;
	t189 = -t190 * t209 + t203 * t216;
	t188 = t192 * t209 - t205 * t216;
	t187 = t192 * t207 + t205 * t215;
	t186 = t190 * t207 + t203 * t215;
	t1 = [t189 * t199 + t193 * t198, -t189 * t198 + t193 * t199, -t186, t189 * t196 + t186 * t201 + t193 * t197 + (t205 * pkin(2) + pkin(8) * t217) * t210 + (-pkin(2) * t217 + t205 * pkin(8)) * t208 + t203 * t218 + t205 * pkin(1) + 0; t188 * t199 + t191 * t198, -t188 * t198 + t191 * t199, t187, t188 * t196 - t187 * t201 + t191 * t197 + (t203 * pkin(2) - pkin(8) * t213) * t210 + (pkin(2) * t213 + t203 * pkin(8)) * t208 - t205 * t218 + t203 * pkin(1) + 0; t195 * t199 - t198 * t214, -t195 * t198 - t199 * t214, t194, t206 * pkin(7) - t194 * t201 + t195 * t196 + qJ(1) + 0 + (pkin(2) * t208 + (-pkin(8) - t197) * t210) * t204; 0, 0, 0, 1;];
	Tc_mdh = t1;
end