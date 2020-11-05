% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPR13 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:41
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRPR13_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR13_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:26
	% EndTime: 2020-11-04 22:41:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:26
	% EndTime: 2020-11-04 22:41:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t123 = cos(qJ(1));
	t122 = sin(qJ(1));
	t1 = [t123, -t122, 0, 0; t122, t123, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:26
	% EndTime: 2020-11-04 22:41:26
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t124 = sin(pkin(6));
	t127 = sin(qJ(1));
	t135 = t127 * t124;
	t126 = sin(qJ(2));
	t134 = t127 * t126;
	t128 = cos(qJ(2));
	t133 = t127 * t128;
	t129 = cos(qJ(1));
	t132 = t129 * t124;
	t131 = t129 * t126;
	t130 = t129 * t128;
	t125 = cos(pkin(6));
	t1 = [-t125 * t134 + t130, -t125 * t133 - t131, t135, t129 * pkin(1) + pkin(8) * t135 + 0; t125 * t131 + t133, t125 * t130 - t134, -t132, t127 * pkin(1) - pkin(8) * t132 + 0; t124 * t126, t124 * t128, t125, t125 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:26
	% EndTime: 2020-11-04 22:41:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t139 = sin(pkin(6));
	t144 = cos(qJ(3));
	t154 = t139 * t144;
	t141 = sin(qJ(3));
	t153 = t141 * t139;
	t142 = sin(qJ(2));
	t152 = t142 * t144;
	t143 = sin(qJ(1));
	t151 = t143 * t142;
	t145 = cos(qJ(2));
	t150 = t143 * t145;
	t146 = cos(qJ(1));
	t149 = t146 * t142;
	t148 = t146 * t145;
	t147 = pkin(2) * t142 - pkin(9) * t145;
	t140 = cos(pkin(6));
	t138 = t145 * pkin(2) + t142 * pkin(9) + pkin(1);
	t137 = -t140 * t149 - t150;
	t136 = t139 * pkin(8) - t147 * t140;
	t1 = [(-t140 * t152 + t153) * t143 + t144 * t148, (t140 * t151 - t148) * t141 + t143 * t154, t140 * t150 + t149, t136 * t143 + t138 * t146 + 0; -t137 * t144 - t146 * t153, t137 * t141 - t146 * t154, -t140 * t148 + t151, -t136 * t146 + t138 * t143 + 0; t139 * t152 + t140 * t141, t140 * t144 - t142 * t153, -t139 * t145, t140 * pkin(8) + t147 * t139 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:26
	% EndTime: 2020-11-04 22:41:26
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (57->34), mult. (112->60), div. (0->0), fcn. (146->10), ass. (0->29)
	t162 = sin(pkin(6));
	t169 = cos(qJ(3));
	t183 = t162 * t169;
	t164 = sin(qJ(4));
	t170 = cos(qJ(2));
	t182 = t164 * t170;
	t165 = sin(qJ(3));
	t181 = t165 * t162;
	t166 = sin(qJ(2));
	t180 = t166 * t169;
	t167 = sin(qJ(1));
	t179 = t167 * t170;
	t168 = cos(qJ(4));
	t178 = t168 * t170;
	t177 = t169 * t170;
	t171 = cos(qJ(1));
	t176 = t171 * t166;
	t175 = t171 * t170;
	t160 = t169 * pkin(3) + t165 * pkin(10) + pkin(2);
	t174 = t170 * pkin(9) - t160 * t166;
	t173 = t165 * pkin(3) - t169 * pkin(10) + pkin(8);
	t163 = cos(pkin(6));
	t158 = t163 * t180 - t181;
	t172 = t158 * t164 + t163 * t178;
	t159 = t164 * t177 - t168 * t166;
	t157 = t162 * t180 + t163 * t165;
	t156 = t166 * pkin(9) + t160 * t170 + pkin(1);
	t155 = t162 * t173 + t174 * t163;
	t1 = [(-t158 * t167 + t169 * t175) * t168 + (t163 * t179 + t176) * t164, -t171 * t159 + t172 * t167, (-t167 * t163 * t166 + t175) * t165 - t167 * t183, t155 * t167 + t156 * t171 + 0; (t158 * t168 - t163 * t182) * t171 + t167 * (t164 * t166 + t168 * t177), -t167 * t159 - t172 * t171, (t163 * t176 + t179) * t165 + t171 * t183, -t155 * t171 + t156 * t167 + 0; t157 * t168 - t162 * t182, -t157 * t164 - t162 * t178, -t163 * t169 + t166 * t181, -t174 * t162 + t173 * t163 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:26
	% EndTime: 2020-11-04 22:41:26
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (83->38), mult. (140->64), div. (0->0), fcn. (174->10), ass. (0->31)
	t196 = sin(qJ(4));
	t200 = cos(qJ(4));
	t190 = pkin(4) * t200 + qJ(5) * t196 + pkin(3);
	t197 = sin(qJ(3));
	t201 = cos(qJ(3));
	t217 = -t201 * pkin(10) + t190 * t197 + pkin(8);
	t194 = sin(pkin(6));
	t215 = t194 * t201;
	t202 = cos(qJ(2));
	t214 = t196 * t202;
	t213 = t197 * t194;
	t198 = sin(qJ(2));
	t212 = t198 * t201;
	t199 = sin(qJ(1));
	t211 = t199 * t202;
	t210 = t200 * t202;
	t209 = t201 * t202;
	t203 = cos(qJ(1));
	t208 = t203 * t198;
	t207 = t203 * t202;
	t186 = t197 * pkin(10) + t190 * t201 + pkin(2);
	t191 = -t196 * pkin(4) + qJ(5) * t200 - pkin(9);
	t205 = t186 * t198 + t191 * t202;
	t195 = cos(pkin(6));
	t188 = t195 * t212 - t213;
	t204 = t188 * t196 + t195 * t210;
	t189 = t196 * t209 - t200 * t198;
	t187 = t194 * t212 + t195 * t197;
	t185 = t186 * t202 - t191 * t198 + pkin(1);
	t184 = -t194 * t217 + t205 * t195;
	t1 = [(-t188 * t199 + t201 * t207) * t200 + (t195 * t211 + t208) * t196, (-t199 * t195 * t198 + t207) * t197 - t199 * t215, t203 * t189 - t204 * t199, -t184 * t199 + t185 * t203 + 0; (t188 * t200 - t195 * t214) * t203 + t199 * (t196 * t198 + t200 * t209), (t195 * t208 + t211) * t197 + t203 * t215, t199 * t189 + t204 * t203, t184 * t203 + t185 * t199 + 0; t187 * t200 - t194 * t214, -t195 * t201 + t198 * t213, t187 * t196 + t194 * t210, t205 * t194 + t217 * t195 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:26
	% EndTime: 2020-11-04 22:41:26
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (130->44), mult. (210->72), div. (0->0), fcn. (270->12), ass. (0->41)
	t239 = sin(qJ(2));
	t244 = cos(qJ(2));
	t237 = sin(qJ(4));
	t242 = cos(qJ(4));
	t247 = pkin(4) + pkin(5);
	t250 = qJ(5) * t242 - t247 * t237 - pkin(9);
	t231 = qJ(5) * t237 + t247 * t242 + pkin(3);
	t238 = sin(qJ(3));
	t243 = cos(qJ(3));
	t246 = pkin(11) - pkin(10);
	t264 = -t231 * t243 + t246 * t238 - pkin(2);
	t267 = -t264 * t239 + t250 * t244;
	t265 = t231 * t238 + t246 * t243 + pkin(8);
	t234 = sin(pkin(6));
	t260 = t234 * t243;
	t235 = cos(pkin(6));
	t259 = t235 * t239;
	t236 = sin(qJ(6));
	t258 = t236 * t244;
	t257 = t238 * t234;
	t256 = t239 * t243;
	t241 = cos(qJ(6));
	t255 = t241 * t244;
	t254 = t243 * t244;
	t228 = t235 * t256 - t257;
	t223 = -t236 * t228 + t235 * t255;
	t225 = t228 * t241 + t235 * t258;
	t249 = t223 * t242 + t237 * t225;
	t248 = -t237 * t223 + t225 * t242;
	t245 = cos(qJ(1));
	t240 = sin(qJ(1));
	t230 = -t239 * t236 + t241 * t254;
	t229 = t236 * t254 + t239 * t241;
	t227 = t234 * t256 + t235 * t238;
	t224 = t227 * t241 + t234 * t258;
	t222 = -t236 * t227 + t234 * t255;
	t221 = t237 * t229 + t230 * t242;
	t220 = t229 * t242 - t237 * t230;
	t219 = -t250 * t239 - t264 * t244 + pkin(1);
	t218 = -t265 * t234 + t267 * t235;
	t1 = [t221 * t245 - t248 * t240, -t245 * t220 - t249 * t240, (t240 * t259 - t245 * t244) * t238 + t240 * t260, -t218 * t240 + t219 * t245 + 0; t221 * t240 + t248 * t245, -t220 * t240 + t249 * t245, (-t240 * t244 - t245 * t259) * t238 - t245 * t260, t218 * t245 + t219 * t240 + 0; -t222 * t237 + t224 * t242, t222 * t242 + t224 * t237, t235 * t243 - t239 * t257, t267 * t234 + t265 * t235 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end