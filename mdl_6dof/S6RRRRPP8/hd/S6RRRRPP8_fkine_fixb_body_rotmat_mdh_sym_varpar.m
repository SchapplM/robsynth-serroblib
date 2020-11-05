% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPP8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:37
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRPP8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:37:10
	% EndTime: 2020-11-04 22:37:10
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:37:10
	% EndTime: 2020-11-04 22:37:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t105 = cos(qJ(1));
	t104 = sin(qJ(1));
	t1 = [t105, -t104, 0, 0; t104, t105, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:37:10
	% EndTime: 2020-11-04 22:37:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t106 = sin(pkin(6));
	t109 = sin(qJ(1));
	t117 = t109 * t106;
	t108 = sin(qJ(2));
	t116 = t109 * t108;
	t110 = cos(qJ(2));
	t115 = t109 * t110;
	t111 = cos(qJ(1));
	t114 = t111 * t106;
	t113 = t111 * t108;
	t112 = t111 * t110;
	t107 = cos(pkin(6));
	t1 = [-t107 * t116 + t112, -t107 * t115 - t113, t117, t111 * pkin(1) + pkin(8) * t117 + 0; t107 * t113 + t115, t107 * t112 - t116, -t114, t109 * pkin(1) - pkin(8) * t114 + 0; t106 * t108, t106 * t110, t107, t107 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:37:10
	% EndTime: 2020-11-04 22:37:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t121 = sin(pkin(6));
	t126 = cos(qJ(3));
	t136 = t121 * t126;
	t123 = sin(qJ(3));
	t135 = t123 * t121;
	t124 = sin(qJ(2));
	t134 = t124 * t126;
	t125 = sin(qJ(1));
	t133 = t125 * t124;
	t127 = cos(qJ(2));
	t132 = t125 * t127;
	t128 = cos(qJ(1));
	t131 = t128 * t124;
	t130 = t128 * t127;
	t129 = pkin(2) * t124 - pkin(9) * t127;
	t122 = cos(pkin(6));
	t120 = t127 * pkin(2) + t124 * pkin(9) + pkin(1);
	t119 = t122 * t131 + t132;
	t118 = t121 * pkin(8) - t129 * t122;
	t1 = [(-t122 * t134 + t135) * t125 + t126 * t130, (t122 * t133 - t130) * t123 + t125 * t136, t122 * t132 + t131, t118 * t125 + t120 * t128 + 0; t119 * t126 - t128 * t135, -t119 * t123 - t128 * t136, -t122 * t130 + t133, -t118 * t128 + t120 * t125 + 0; t121 * t134 + t122 * t123, t122 * t126 - t124 * t135, -t121 * t127, t122 * pkin(8) + t129 * t121 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:37:10
	% EndTime: 2020-11-04 22:37:10
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (57->34), mult. (112->60), div. (0->0), fcn. (146->10), ass. (0->29)
	t144 = sin(pkin(6));
	t151 = cos(qJ(3));
	t165 = t144 * t151;
	t146 = sin(qJ(4));
	t152 = cos(qJ(2));
	t164 = t146 * t152;
	t147 = sin(qJ(3));
	t163 = t147 * t144;
	t148 = sin(qJ(2));
	t162 = t148 * t151;
	t149 = sin(qJ(1));
	t161 = t149 * t152;
	t150 = cos(qJ(4));
	t160 = t150 * t152;
	t159 = t151 * t152;
	t153 = cos(qJ(1));
	t158 = t153 * t148;
	t157 = t153 * t152;
	t142 = t151 * pkin(3) + t147 * pkin(10) + pkin(2);
	t156 = t152 * pkin(9) - t142 * t148;
	t155 = t147 * pkin(3) - t151 * pkin(10) + pkin(8);
	t145 = cos(pkin(6));
	t140 = t145 * t162 - t163;
	t154 = t140 * t146 + t145 * t160;
	t141 = t146 * t159 - t150 * t148;
	t139 = t144 * t162 + t145 * t147;
	t138 = t148 * pkin(9) + t142 * t152 + pkin(1);
	t137 = t144 * t155 + t156 * t145;
	t1 = [(-t140 * t149 + t151 * t157) * t150 + (t145 * t161 + t158) * t146, -t153 * t141 + t154 * t149, (-t149 * t145 * t148 + t157) * t147 - t149 * t165, t137 * t149 + t138 * t153 + 0; (t140 * t150 - t145 * t164) * t153 + t149 * (t148 * t146 + t150 * t159), -t149 * t141 - t154 * t153, (t145 * t158 + t161) * t147 + t153 * t165, -t137 * t153 + t138 * t149 + 0; t139 * t150 - t144 * t164, -t139 * t146 - t144 * t160, -t145 * t151 + t148 * t163, -t156 * t144 + t155 * t145 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:37:10
	% EndTime: 2020-11-04 22:37:10
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (83->38), mult. (140->64), div. (0->0), fcn. (174->10), ass. (0->31)
	t178 = sin(qJ(4));
	t182 = cos(qJ(4));
	t172 = pkin(4) * t182 + qJ(5) * t178 + pkin(3);
	t179 = sin(qJ(3));
	t183 = cos(qJ(3));
	t199 = -t183 * pkin(10) + t172 * t179 + pkin(8);
	t176 = sin(pkin(6));
	t197 = t176 * t183;
	t184 = cos(qJ(2));
	t196 = t178 * t184;
	t195 = t179 * t176;
	t180 = sin(qJ(2));
	t194 = t180 * t183;
	t181 = sin(qJ(1));
	t193 = t181 * t184;
	t192 = t182 * t184;
	t191 = t183 * t184;
	t185 = cos(qJ(1));
	t190 = t185 * t180;
	t189 = t185 * t184;
	t168 = t179 * pkin(10) + t172 * t183 + pkin(2);
	t173 = -t178 * pkin(4) + qJ(5) * t182 - pkin(9);
	t187 = t168 * t180 + t173 * t184;
	t177 = cos(pkin(6));
	t170 = t177 * t194 - t195;
	t186 = t170 * t178 + t177 * t192;
	t171 = t178 * t191 - t182 * t180;
	t169 = t176 * t194 + t177 * t179;
	t167 = t168 * t184 - t173 * t180 + pkin(1);
	t166 = -t199 * t176 + t187 * t177;
	t1 = [(-t170 * t181 + t183 * t189) * t182 + (t177 * t193 + t190) * t178, (-t181 * t177 * t180 + t189) * t179 - t181 * t197, t185 * t171 - t186 * t181, -t166 * t181 + t167 * t185 + 0; (t170 * t182 - t177 * t196) * t185 + t181 * (t180 * t178 + t182 * t191), (t177 * t190 + t193) * t179 + t185 * t197, t181 * t171 + t186 * t185, t166 * t185 + t167 * t181 + 0; t169 * t182 - t176 * t196, -t177 * t183 + t180 * t195, t169 * t178 + t176 * t192, t187 * t176 + t199 * t177 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:37:10
	% EndTime: 2020-11-04 22:37:10
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (104->40), mult. (142->64), div. (0->0), fcn. (176->10), ass. (0->33)
	t213 = sin(qJ(4));
	t217 = cos(qJ(4));
	t221 = pkin(4) + pkin(5);
	t206 = qJ(5) * t213 + t221 * t217 + pkin(3);
	t212 = qJ(6) - pkin(10);
	t214 = sin(qJ(3));
	t218 = cos(qJ(3));
	t202 = t206 * t218 - t212 * t214 + pkin(2);
	t215 = sin(qJ(2));
	t219 = cos(qJ(2));
	t223 = qJ(5) * t217 - t221 * t213 - pkin(9);
	t239 = t202 * t215 + t223 * t219;
	t237 = t206 * t214 + t212 * t218 + pkin(8);
	t210 = sin(pkin(6));
	t234 = t210 * t218;
	t233 = t213 * t219;
	t232 = t214 * t210;
	t231 = t215 * t218;
	t216 = sin(qJ(1));
	t230 = t216 * t219;
	t229 = t217 * t219;
	t228 = t218 * t219;
	t220 = cos(qJ(1));
	t227 = t220 * t215;
	t226 = t220 * t219;
	t211 = cos(pkin(6));
	t204 = t211 * t231 - t232;
	t222 = t204 * t213 + t211 * t229;
	t205 = t213 * t228 - t217 * t215;
	t203 = t210 * t231 + t211 * t214;
	t201 = t202 * t219 - t223 * t215 + pkin(1);
	t200 = -t237 * t210 + t239 * t211;
	t1 = [(-t204 * t216 + t218 * t226) * t217 + (t211 * t230 + t227) * t213, t220 * t205 - t222 * t216, (t216 * t211 * t215 - t226) * t214 + t216 * t234, -t200 * t216 + t201 * t220 + 0; (t204 * t217 - t211 * t233) * t220 + t216 * (t215 * t213 + t217 * t228), t216 * t205 + t222 * t220, (-t211 * t227 - t230) * t214 - t220 * t234, t200 * t220 + t201 * t216 + 0; t203 * t217 - t210 * t233, t203 * t213 + t210 * t229, t211 * t218 - t215 * t232, t239 * t210 + t237 * t211 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end