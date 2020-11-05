% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:03
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:03:39
	% EndTime: 2020-11-04 21:03:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:03:39
	% EndTime: 2020-11-04 21:03:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t113 = cos(pkin(11));
	t112 = sin(pkin(11));
	t1 = [t113, -t112, 0, 0; t112, t113, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:03:39
	% EndTime: 2020-11-04 21:03:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t114 = sin(pkin(11));
	t115 = sin(pkin(6));
	t123 = t114 * t115;
	t116 = cos(pkin(11));
	t122 = t116 * t115;
	t117 = cos(pkin(6));
	t118 = sin(qJ(2));
	t121 = t117 * t118;
	t119 = cos(qJ(2));
	t120 = t117 * t119;
	t1 = [-t114 * t121 + t116 * t119, -t114 * t120 - t116 * t118, t123, t116 * pkin(1) + pkin(7) * t123 + 0; t114 * t119 + t116 * t121, -t114 * t118 + t116 * t120, -t122, t114 * pkin(1) - pkin(7) * t122 + 0; t115 * t118, t115 * t119, t117, t117 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:03:39
	% EndTime: 2020-11-04 21:03:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (53->22), mult. (37->26), div. (0->0), fcn. (44->12), ass. (0->20)
	t143 = pkin(2) * sin(qJ(2));
	t136 = qJ(2) + pkin(12);
	t141 = pkin(7) + qJ(3);
	t140 = cos(pkin(6));
	t139 = cos(pkin(11));
	t138 = sin(pkin(6));
	t137 = sin(pkin(11));
	t135 = cos(t136);
	t134 = sin(t136);
	t133 = pkin(6) - t136;
	t132 = pkin(6) + t136;
	t131 = cos(qJ(2)) * pkin(2) + pkin(1);
	t130 = cos(t132);
	t129 = sin(t133);
	t128 = cos(t133) / 0.2e1;
	t127 = sin(t132) / 0.2e1;
	t126 = -t138 * t141 + t140 * t143;
	t125 = t130 / 0.2e1 + t128;
	t124 = t127 - t129 / 0.2e1;
	t1 = [-t137 * t124 + t139 * t135, -t137 * t125 - t139 * t134, t137 * t138, -t137 * t126 + t139 * t131 + 0; t139 * t124 + t137 * t135, t139 * t125 - t137 * t134, -t139 * t138, t139 * t126 + t137 * t131 + 0; t128 - t130 / 0.2e1, t129 / 0.2e1 + t127, t140, t138 * t143 + t140 * t141 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:03:39
	% EndTime: 2020-11-04 21:03:39
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (73->35), mult. (85->53), div. (0->0), fcn. (105->16), ass. (0->27)
	t155 = sin(pkin(11));
	t159 = cos(pkin(6));
	t169 = t155 * t159;
	t156 = sin(pkin(6));
	t160 = pkin(7) + qJ(3);
	t168 = t156 * t160;
	t161 = sin(qJ(4));
	t167 = t156 * t161;
	t163 = cos(qJ(4));
	t166 = t156 * t163;
	t158 = cos(pkin(11));
	t165 = t158 * t159;
	t153 = qJ(2) + pkin(12);
	t164 = cos(qJ(2));
	t162 = sin(qJ(2));
	t157 = cos(pkin(12));
	t154 = sin(pkin(12));
	t152 = cos(t153);
	t151 = sin(t153);
	t150 = pkin(6) - t153;
	t149 = pkin(6) + t153;
	t148 = -t154 * pkin(3) + t157 * pkin(8);
	t147 = t157 * pkin(3) + t154 * pkin(8) + pkin(2);
	t146 = cos(t149) + cos(t150);
	t145 = t151 * t165 + t155 * t152;
	t144 = t151 * t169 - t158 * t152;
	t1 = [-t144 * t163 + t155 * t167, t144 * t161 + t155 * t166, t158 * t151 + t155 * t146 / 0.2e1, (t158 * t147 + t148 * t169) * t164 + (-t147 * t169 + t158 * t148) * t162 + t155 * t168 + t158 * pkin(1) + 0; t145 * t163 - t158 * t167, -t145 * t161 - t158 * t166, t155 * t151 - t158 * t146 / 0.2e1, (t155 * t147 - t148 * t165) * t164 + (t147 * t165 + t155 * t148) * t162 - t158 * t168 + t155 * pkin(1) + 0; t151 * t166 + t159 * t161, -t151 * t167 + t159 * t163, -sin(t150) / 0.2e1 - sin(t149) / 0.2e1, t159 * t160 + qJ(1) + 0 + (t147 * t162 - t148 * t164) * t156; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:03:39
	% EndTime: 2020-11-04 21:03:39
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (108->49), mult. (117->69), div. (0->0), fcn. (141->18), ass. (0->31)
	t204 = pkin(4) * cos(qJ(4));
	t184 = qJ(2) + pkin(12);
	t179 = sin(t184);
	t188 = sin(pkin(6));
	t203 = t179 * t188;
	t187 = sin(pkin(11));
	t202 = t187 * t188;
	t190 = cos(pkin(11));
	t201 = t188 * t190;
	t191 = cos(pkin(6));
	t200 = t191 * t187;
	t199 = t191 * t190;
	t186 = sin(pkin(12));
	t189 = cos(pkin(12));
	t197 = pkin(9) + pkin(8);
	t173 = t189 * pkin(3) + t197 * t186 + pkin(2);
	t198 = pkin(4) * sin(qJ(4)) + pkin(7) + qJ(3);
	t196 = cos(qJ(2));
	t194 = sin(qJ(2));
	t185 = qJ(4) + qJ(5);
	t182 = cos(t185);
	t181 = sin(t185);
	t180 = cos(t184);
	t178 = pkin(6) - t184;
	t177 = pkin(6) + t184;
	t176 = t197 * t189;
	t174 = -t186 * pkin(3) + t176;
	t172 = cos(t177) + cos(t178);
	t171 = t179 * t199 + t187 * t180;
	t170 = t179 * t200 - t190 * t180;
	t1 = [-t170 * t182 + t181 * t202, t170 * t181 + t182 * t202, t190 * t179 + t187 * t172 / 0.2e1, (-(t186 * t200 - t190 * t189) * t204 + t174 * t200 + t173 * t190) * t196 + (-(t190 * t186 + t189 * t200) * t204 - t173 * t200 + t174 * t190) * t194 + t190 * pkin(1) + 0 + t198 * t202; t171 * t182 - t181 * t201, -t171 * t181 - t182 * t201, t187 * t179 - t190 * t172 / 0.2e1, ((t186 * t199 + t189 * t187) * t204 - t174 * t199 + t173 * t187) * t196 + ((-t187 * t186 + t189 * t199) * t204 + t173 * t199 + t174 * t187) * t194 + t187 * pkin(1) + 0 - t198 * t201; t191 * t181 + t182 * t203, -t181 * t203 + t191 * t182, -sin(t178) / 0.2e1 - sin(t177) / 0.2e1, qJ(1) + 0 + t198 * t191 + ((t189 * t204 + t173) * t194 - (t176 + (-pkin(3) - t204) * t186) * t196) * t188; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:03:39
	% EndTime: 2020-11-04 21:03:39
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (187->62), mult. (188->87), div. (0->0), fcn. (234->20), ass. (0->42)
	t250 = pkin(4) * cos(qJ(4));
	t228 = qJ(2) + pkin(12);
	t223 = sin(t228);
	t232 = sin(pkin(6));
	t249 = t223 * t232;
	t231 = sin(pkin(11));
	t248 = t231 * t232;
	t234 = cos(pkin(11));
	t247 = t232 * t234;
	t235 = cos(pkin(6));
	t246 = t234 * t235;
	t245 = t235 * t231;
	t230 = sin(pkin(12));
	t233 = cos(pkin(12));
	t243 = pkin(9) + pkin(8);
	t217 = t233 * pkin(3) + t243 * t230 + pkin(2);
	t244 = pkin(4) * sin(qJ(4)) + pkin(7) + qJ(3);
	t242 = cos(qJ(2));
	t240 = cos(qJ(6));
	t239 = sin(qJ(2));
	t237 = sin(qJ(6));
	t229 = qJ(4) + qJ(5);
	t226 = cos(t229);
	t225 = sin(t229);
	t224 = cos(t228);
	t222 = pkin(6) - t228;
	t221 = pkin(6) + t228;
	t220 = t243 * t233;
	t218 = -t230 * pkin(3) + t220;
	t216 = cos(t221) + cos(t222);
	t215 = -sin(t222) / 0.2e1 - sin(t221) / 0.2e1;
	t214 = t223 * t246 + t231 * t224;
	t213 = t223 * t245 - t234 * t224;
	t212 = t235 * t225 + t226 * t249;
	t211 = t225 * t249 - t235 * t226;
	t210 = t234 * t223 + t231 * t216 / 0.2e1;
	t209 = t231 * t223 - t234 * t216 / 0.2e1;
	t208 = -t213 * t226 + t225 * t248;
	t207 = t214 * t226 - t225 * t247;
	t206 = t214 * t225 + t226 * t247;
	t205 = t213 * t225 + t226 * t248;
	t1 = [t208 * t240 + t210 * t237, -t208 * t237 + t210 * t240, -t205, t208 * pkin(5) - t205 * pkin(10) + (-(t230 * t245 - t234 * t233) * t250 + t218 * t245 + t217 * t234) * t242 + (-(t234 * t230 + t233 * t245) * t250 - t217 * t245 + t218 * t234) * t239 + t234 * pkin(1) + 0 + t244 * t248; t207 * t240 + t209 * t237, -t207 * t237 + t209 * t240, t206, t207 * pkin(5) + t206 * pkin(10) + ((t230 * t246 + t233 * t231) * t250 - t218 * t246 + t217 * t231) * t242 + ((-t231 * t230 + t233 * t246) * t250 + t217 * t246 + t218 * t231) * t239 + t231 * pkin(1) + 0 - t244 * t247; t212 * t240 + t215 * t237, -t212 * t237 + t215 * t240, t211, t212 * pkin(5) + t211 * pkin(10) + qJ(1) + 0 + t244 * t235 + ((t233 * t250 + t217) * t239 - (t220 + (-pkin(3) - t250) * t230) * t242) * t232; 0, 0, 0, 1;];
	Tc_mdh = t1;
end