% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRR8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:32
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:32:10
	% EndTime: 2020-11-04 22:32:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:32:10
	% EndTime: 2020-11-04 22:32:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t99 = cos(qJ(1));
	t98 = sin(qJ(1));
	t1 = [t99, -t98, 0, 0; t98, t99, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:32:10
	% EndTime: 2020-11-04 22:32:10
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t100 = sin(pkin(6));
	t103 = sin(qJ(1));
	t111 = t103 * t100;
	t102 = sin(qJ(2));
	t110 = t103 * t102;
	t104 = cos(qJ(2));
	t109 = t103 * t104;
	t105 = cos(qJ(1));
	t108 = t105 * t100;
	t107 = t105 * t102;
	t106 = t105 * t104;
	t101 = cos(pkin(6));
	t1 = [-t101 * t110 + t106, -t101 * t109 - t107, t111, t105 * pkin(1) + pkin(8) * t111 + 0; t101 * t107 + t109, t101 * t106 - t110, -t108, t103 * pkin(1) - pkin(8) * t108 + 0; t100 * t102, t100 * t104, t101, t101 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:32:10
	% EndTime: 2020-11-04 22:32:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t115 = sin(pkin(6));
	t120 = cos(qJ(3));
	t130 = t115 * t120;
	t117 = sin(qJ(3));
	t129 = t117 * t115;
	t118 = sin(qJ(2));
	t128 = t118 * t120;
	t119 = sin(qJ(1));
	t127 = t119 * t118;
	t121 = cos(qJ(2));
	t126 = t119 * t121;
	t122 = cos(qJ(1));
	t125 = t122 * t118;
	t124 = t122 * t121;
	t123 = pkin(2) * t118 - pkin(9) * t121;
	t116 = cos(pkin(6));
	t114 = t121 * pkin(2) + t118 * pkin(9) + pkin(1);
	t113 = -t116 * t125 - t126;
	t112 = t115 * pkin(8) - t123 * t116;
	t1 = [(-t116 * t128 + t129) * t119 + t120 * t124, (t116 * t127 - t124) * t117 + t119 * t130, t116 * t126 + t125, t112 * t119 + t114 * t122 + 0; -t113 * t120 - t122 * t129, t113 * t117 - t122 * t130, -t116 * t124 + t127, -t112 * t122 + t114 * t119 + 0; t115 * t128 + t116 * t117, t116 * t120 - t118 * t129, -t115 * t121, t116 * pkin(8) + t123 * t115 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:32:10
	% EndTime: 2020-11-04 22:32:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (54->26), mult. (69->38), div. (0->0), fcn. (90->10), ass. (0->25)
	t139 = sin(pkin(6));
	t143 = sin(qJ(2));
	t155 = t139 * t143;
	t144 = sin(qJ(1));
	t154 = t139 * t144;
	t146 = cos(qJ(1));
	t153 = t139 * t146;
	t152 = t144 * t143;
	t145 = cos(qJ(2));
	t151 = t144 * t145;
	t150 = t146 * t143;
	t149 = t146 * t145;
	t148 = sin(qJ(3)) * pkin(3) + pkin(8);
	t135 = cos(qJ(3)) * pkin(3) + pkin(2);
	t141 = qJ(4) + pkin(9);
	t147 = t135 * t143 - t141 * t145;
	t140 = cos(pkin(6));
	t138 = qJ(3) + pkin(12);
	t137 = cos(t138);
	t136 = sin(t138);
	t134 = t140 * t150 + t151;
	t133 = t140 * t152 - t149;
	t132 = t135 * t145 + t141 * t143 + pkin(1);
	t131 = t139 * t148 - t147 * t140;
	t1 = [-t133 * t137 + t136 * t154, t133 * t136 + t137 * t154, t140 * t151 + t150, t131 * t144 + t132 * t146 + 0; t134 * t137 - t136 * t153, -t134 * t136 - t137 * t153, -t140 * t149 + t152, -t131 * t146 + t132 * t144 + 0; t140 * t136 + t137 * t155, -t136 * t155 + t140 * t137, -t139 * t145, t147 * t139 + t148 * t140 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:32:10
	% EndTime: 2020-11-04 22:32:10
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (104->39), mult. (152->61), div. (0->0), fcn. (186->14), ass. (0->40)
	t173 = sin(pkin(12));
	t175 = cos(pkin(12));
	t167 = t175 * pkin(4) + t173 * pkin(10) + pkin(3);
	t179 = sin(qJ(3));
	t164 = t167 * t179 + pkin(8);
	t168 = -t173 * pkin(4) + t175 * pkin(10);
	t165 = t168 * t179 + pkin(2);
	t174 = sin(pkin(6));
	t176 = cos(pkin(6));
	t180 = sin(qJ(2));
	t183 = cos(qJ(3));
	t177 = qJ(4) + pkin(9);
	t184 = cos(qJ(2));
	t194 = t177 * t184;
	t199 = (t165 * t180 - t194) * t176 + (t176 * t167 * t180 + t174 * t168) * t183 - t164 * t174;
	t198 = t174 * t180;
	t181 = sin(qJ(1));
	t197 = t174 * t181;
	t196 = t174 * t184;
	t185 = cos(qJ(1));
	t195 = t174 * t185;
	t193 = t181 * t180;
	t192 = t181 * t184;
	t191 = t185 * t180;
	t190 = t185 * t184;
	t160 = t176 * t193 - t190;
	t172 = qJ(3) + pkin(12);
	t170 = sin(t172);
	t171 = cos(t172);
	t187 = -t160 * t171 + t170 * t197;
	t162 = t176 * t191 + t192;
	t186 = -t162 * t171 + t170 * t195;
	t182 = cos(qJ(5));
	t178 = sin(qJ(5));
	t163 = t176 * t192 + t191;
	t161 = t176 * t190 - t193;
	t159 = t176 * t170 + t171 * t198;
	t158 = t167 * t183 + t165;
	t156 = t158 * t184 + t177 * t180 + pkin(1);
	t1 = [t163 * t178 + t187 * t182, t163 * t182 - t187 * t178, -t160 * t170 - t171 * t197, t156 * t185 - t199 * t181 + 0; -t178 * t161 - t186 * t182, -t182 * t161 + t186 * t178, t162 * t170 + t171 * t195, t156 * t181 + t199 * t185 + 0; t159 * t182 - t178 * t196, -t159 * t178 - t182 * t196, t170 * t198 - t176 * t171, pkin(7) + 0 + (t158 * t180 - t194) * t174 + (-t168 * t183 + t164) * t176; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:32:10
	% EndTime: 2020-11-04 22:32:10
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (123->45), mult. (155->61), div. (0->0), fcn. (201->14), ass. (0->40)
	t240 = pkin(5) * sin(qJ(5));
	t213 = cos(qJ(3)) * pkin(3) + pkin(2);
	t225 = sin(qJ(2));
	t239 = t213 * t225;
	t220 = sin(pkin(6));
	t238 = t220 * t225;
	t226 = sin(qJ(1));
	t237 = t220 * t226;
	t227 = cos(qJ(2));
	t236 = t220 * t227;
	t228 = cos(qJ(1));
	t235 = t220 * t228;
	t234 = t226 * t225;
	t233 = t226 * t227;
	t232 = t228 * t225;
	t231 = t228 * t227;
	t230 = sin(qJ(3)) * pkin(3) + pkin(8);
	t229 = -pkin(11) - pkin(10);
	t222 = qJ(4) + pkin(9);
	t221 = cos(pkin(6));
	t219 = qJ(5) + qJ(6);
	t218 = qJ(3) + pkin(12);
	t217 = cos(t219);
	t216 = sin(t219);
	t215 = cos(t218);
	t214 = sin(t218);
	t212 = cos(qJ(5)) * pkin(5) + pkin(4);
	t211 = t221 * t233 + t232;
	t210 = t221 * t232 + t233;
	t209 = -t221 * t231 + t234;
	t208 = t221 * t234 - t231;
	t207 = t213 * t227 + t222 * t225 + pkin(1);
	t206 = t221 * t214 + t215 * t238;
	t205 = t214 * t238 - t221 * t215;
	t204 = t220 * t230 + (t222 * t227 - t239) * t221;
	t203 = -t208 * t215 + t214 * t237;
	t202 = t210 * t215 - t214 * t235;
	t201 = t210 * t214 + t215 * t235;
	t200 = t208 * t214 + t215 * t237;
	t1 = [t203 * t217 + t211 * t216, -t203 * t216 + t211 * t217, -t200, t200 * t229 + t203 * t212 + t204 * t226 + t207 * t228 + t211 * t240 + 0; t202 * t217 + t209 * t216, -t202 * t216 + t209 * t217, t201, -t201 * t229 + t202 * t212 - t204 * t228 + t207 * t226 + t209 * t240 + 0; t206 * t217 - t216 * t236, -t206 * t216 - t217 * t236, t205, -t205 * t229 + t206 * t212 + pkin(7) + 0 + t230 * t221 + (t239 + (-t222 - t240) * t227) * t220; 0, 0, 0, 1;];
	Tc_mdh = t1;
end