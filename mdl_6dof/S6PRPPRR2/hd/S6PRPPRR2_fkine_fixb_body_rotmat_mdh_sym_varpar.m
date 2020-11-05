% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPPRR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:57
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:57:44
	% EndTime: 2020-11-04 20:57:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:57:44
	% EndTime: 2020-11-04 20:57:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t128 = cos(pkin(10));
	t127 = sin(pkin(10));
	t1 = [t128, -t127, 0, 0; t127, t128, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:57:44
	% EndTime: 2020-11-04 20:57:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t129 = sin(pkin(10));
	t130 = sin(pkin(6));
	t138 = t129 * t130;
	t131 = cos(pkin(10));
	t137 = t131 * t130;
	t132 = cos(pkin(6));
	t133 = sin(qJ(2));
	t136 = t132 * t133;
	t134 = cos(qJ(2));
	t135 = t132 * t134;
	t1 = [-t129 * t136 + t131 * t134, -t129 * t135 - t131 * t133, t138, t131 * pkin(1) + pkin(7) * t138 + 0; t129 * t134 + t131 * t136, -t129 * t133 + t131 * t135, -t137, t129 * pkin(1) - pkin(7) * t137 + 0; t130 * t133, t130 * t134, t132, t132 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:57:44
	% EndTime: 2020-11-04 20:57:44
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (53->22), mult. (37->26), div. (0->0), fcn. (44->12), ass. (0->20)
	t158 = pkin(2) * sin(qJ(2));
	t151 = qJ(2) + pkin(11);
	t156 = pkin(7) + qJ(3);
	t155 = cos(pkin(6));
	t154 = cos(pkin(10));
	t153 = sin(pkin(6));
	t152 = sin(pkin(10));
	t150 = cos(t151);
	t149 = sin(t151);
	t148 = pkin(6) - t151;
	t147 = pkin(6) + t151;
	t146 = cos(qJ(2)) * pkin(2) + pkin(1);
	t145 = cos(t147);
	t144 = sin(t148);
	t143 = cos(t148) / 0.2e1;
	t142 = sin(t147) / 0.2e1;
	t141 = -t153 * t156 + t155 * t158;
	t140 = t145 / 0.2e1 + t143;
	t139 = t142 - t144 / 0.2e1;
	t1 = [-t152 * t139 + t154 * t150, -t152 * t140 - t154 * t149, t152 * t153, -t152 * t141 + t154 * t146 + 0; t154 * t139 + t152 * t150, t154 * t140 - t152 * t149, -t154 * t153, t154 * t141 + t152 * t146 + 0; t143 - t145 / 0.2e1, t144 / 0.2e1 + t142, t155, t153 * t158 + t155 * t156 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:57:44
	% EndTime: 2020-11-04 20:57:44
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (73->32), mult. (65->44), div. (0->0), fcn. (76->14), ass. (0->27)
	t173 = sin(pkin(10));
	t174 = sin(pkin(6));
	t184 = t173 * t174;
	t177 = cos(pkin(6));
	t183 = t173 * t177;
	t176 = cos(pkin(10));
	t182 = t176 * t174;
	t181 = t176 * t177;
	t171 = qJ(2) + pkin(11);
	t180 = cos(qJ(2));
	t179 = sin(qJ(2));
	t178 = pkin(7) + qJ(3);
	t175 = cos(pkin(11));
	t172 = sin(pkin(11));
	t170 = cos(t171);
	t169 = sin(t171);
	t168 = pkin(6) - t171;
	t167 = pkin(6) + t171;
	t166 = cos(t168);
	t165 = cos(t167);
	t164 = sin(t168);
	t163 = sin(t167);
	t162 = -t172 * pkin(3) + qJ(4) * t175;
	t161 = t175 * pkin(3) + qJ(4) * t172 + pkin(2);
	t160 = t165 + t166;
	t159 = -t163 + t164;
	t1 = [t184, -t176 * t170 - t173 * t159 / 0.2e1, t176 * t169 + t173 * t160 / 0.2e1, (t176 * t161 + t162 * t183) * t180 + (-t161 * t183 + t176 * t162) * t179 + t178 * t184 + t176 * pkin(1) + 0; -t182, -t173 * t170 + t176 * t159 / 0.2e1, t173 * t169 - t176 * t160 / 0.2e1, (t173 * t161 - t162 * t181) * t180 + (t161 * t181 + t173 * t162) * t179 - t178 * t182 + t173 * pkin(1) + 0; t177, -t166 / 0.2e1 + t165 / 0.2e1, -t164 / 0.2e1 - t163 / 0.2e1, t177 * t178 + qJ(1) + 0 + (t161 * t179 - t162 * t180) * t174; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:57:44
	% EndTime: 2020-11-04 20:57:44
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (98->41), mult. (121->60), div. (0->0), fcn. (119->16), ass. (0->33)
	t199 = qJ(2) + pkin(11);
	t197 = cos(t199);
	t201 = sin(pkin(10));
	t218 = t201 * t197;
	t205 = cos(pkin(6));
	t217 = t201 * t205;
	t198 = qJ(3) + pkin(4) + pkin(7);
	t202 = sin(pkin(6));
	t216 = t202 * t198;
	t206 = sin(qJ(5));
	t215 = t202 * t206;
	t208 = cos(qJ(5));
	t214 = t202 * t208;
	t204 = cos(pkin(10));
	t213 = t204 * t197;
	t212 = t204 * t205;
	t194 = pkin(6) + t199;
	t190 = sin(t194);
	t195 = pkin(6) - t199;
	t191 = sin(t195);
	t211 = t191 / 0.2e1 - t190 / 0.2e1;
	t210 = cos(t195) / 0.2e1 - cos(t194) / 0.2e1;
	t209 = cos(qJ(2));
	t207 = sin(qJ(2));
	t203 = cos(pkin(11));
	t200 = sin(pkin(11));
	t196 = sin(t199);
	t189 = -t200 * pkin(3) + qJ(4) * t203;
	t188 = t203 * pkin(3) + qJ(4) * t200 + pkin(2);
	t187 = -t190 + t191;
	t186 = t204 * t196 + t197 * t217;
	t185 = t201 * t196 - t197 * t212;
	t1 = [t186 * t206 + t201 * t214, t186 * t208 - t201 * t215, t213 + t201 * t187 / 0.2e1, t201 * t216 + t204 * pkin(1) + 0 + (t204 * t188 + t189 * t217) * t209 + (-t188 * t217 + t204 * t189) * t207 + (t211 * t201 + t213) * pkin(8); t185 * t206 - t204 * t214, t185 * t208 + t204 * t215, t218 - t204 * t187 / 0.2e1, -t204 * t216 + t201 * pkin(1) + 0 + (t201 * t188 - t189 * t212) * t209 + (t188 * t212 + t201 * t189) * t207 + (-t211 * t204 + t218) * pkin(8); -t197 * t215 + t205 * t208, -t197 * t214 - t205 * t206, t210, t198 * t205 + 0 + qJ(1) + (t188 * t207 - t189 * t209) * t202 + t210 * pkin(8); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:57:44
	% EndTime: 2020-11-04 20:57:44
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (118->65), mult. (194->100), div. (0->0), fcn. (230->16), ass. (0->42)
	t234 = sin(pkin(6));
	t239 = sin(qJ(5));
	t242 = cos(qJ(5));
	t245 = pkin(5) * t242 + pkin(9) * t239 + pkin(4) + pkin(7) + qJ(3);
	t268 = t245 * t234;
	t267 = -t239 * pkin(5) + t242 * pkin(9);
	t232 = sin(pkin(11));
	t264 = qJ(4) * t232 + pkin(2);
	t231 = qJ(2) + pkin(11);
	t226 = sin(t231);
	t263 = t226 * t234;
	t233 = sin(pkin(10));
	t237 = cos(pkin(6));
	t262 = t233 * t237;
	t241 = cos(qJ(6));
	t261 = t233 * t241;
	t260 = t234 * t239;
	t259 = t234 * t242;
	t235 = cos(pkin(11));
	t258 = t235 * t233;
	t236 = cos(pkin(10));
	t257 = t236 * t237;
	t238 = sin(qJ(6));
	t256 = t238 * t233;
	t255 = t238 * t236;
	t254 = t238 * t239;
	t253 = t241 * t236;
	t252 = t233 * t254;
	t251 = t239 * t261;
	t250 = t236 * t254;
	t249 = t237 * t253;
	t248 = t233 * t259;
	t247 = t236 * t259;
	t244 = pkin(3) + pkin(8);
	t243 = cos(qJ(2));
	t240 = sin(qJ(2));
	t229 = qJ(4) * t235;
	t227 = cos(t231);
	t225 = -t232 * t244 + t229;
	t224 = t244 * t235 + t264;
	t219 = -t227 * t260 + t237 * t242;
	t1 = [(t237 * t251 + t255) * t227 + (-t237 * t256 + t239 * t253) * t226 + t241 * t248, (-t237 * t252 + t253) * t227 + (-t237 * t261 - t250) * t226 - t238 * t248, t233 * t260 - (t236 * t226 + t227 * t262) * t242, (t224 * t236 - t267 * (t236 * t232 + t237 * t258)) * t243 + (t225 * t236 + t267 * (t232 * t262 - t236 * t235)) * t240 + t236 * pkin(1) + 0 + ((-t224 * t240 + t225 * t243) * t237 + t268) * t233; (-t239 * t249 + t256) * t227 + (t237 * t255 + t251) * t226 - t241 * t247, (t237 * t250 + t261) * t227 + (t249 - t252) * t226 + t238 * t247, -t236 * t260 + (-t233 * t226 + t227 * t257) * t242, (-t225 * t257 + t224 * t233 + t267 * (-t233 * t232 + t235 * t257)) * t243 + (t224 * t257 + t225 * t233 - t267 * (t232 * t257 + t258)) * t240 + t233 * pkin(1) + 0 - t236 * t268; t219 * t241 + t238 * t263, -t219 * t238 + t241 * t263, t227 * t259 + t237 * t239, qJ(1) + 0 + (cos(pkin(6) - t231) / 0.2e1 - cos(pkin(6) + t231) / 0.2e1) * pkin(8) + t245 * t237 + (t267 * t227 + (t235 * pkin(3) + t264) * t240 - (-t232 * pkin(3) + t229) * t243) * t234; 0, 0, 0, 1;];
	Tc_mdh = t1;
end