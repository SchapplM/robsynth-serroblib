% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:48
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:08
	% EndTime: 2020-11-04 22:48:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:08
	% EndTime: 2020-11-04 22:48:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t101 = cos(qJ(1));
	t100 = sin(qJ(1));
	t1 = [t101, -t100, 0, 0; t100, t101, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:08
	% EndTime: 2020-11-04 22:48:08
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t102 = sin(pkin(6));
	t105 = sin(qJ(1));
	t113 = t105 * t102;
	t104 = sin(qJ(2));
	t112 = t105 * t104;
	t106 = cos(qJ(2));
	t111 = t105 * t106;
	t107 = cos(qJ(1));
	t110 = t107 * t102;
	t109 = t107 * t104;
	t108 = t107 * t106;
	t103 = cos(pkin(6));
	t1 = [-t103 * t112 + t108, -t103 * t111 - t109, t113, t107 * pkin(1) + pkin(8) * t113 + 0; t103 * t109 + t111, t103 * t108 - t112, -t110, t105 * pkin(1) - pkin(8) * t110 + 0; t102 * t104, t102 * t106, t103, t103 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:08
	% EndTime: 2020-11-04 22:48:08
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t117 = sin(pkin(6));
	t122 = cos(qJ(3));
	t132 = t117 * t122;
	t119 = sin(qJ(3));
	t131 = t119 * t117;
	t120 = sin(qJ(2));
	t130 = t120 * t122;
	t121 = sin(qJ(1));
	t129 = t121 * t120;
	t123 = cos(qJ(2));
	t128 = t121 * t123;
	t124 = cos(qJ(1));
	t127 = t124 * t120;
	t126 = t124 * t123;
	t125 = pkin(2) * t120 - pkin(9) * t123;
	t118 = cos(pkin(6));
	t116 = t123 * pkin(2) + t120 * pkin(9) + pkin(1);
	t115 = -t118 * t127 - t128;
	t114 = t117 * pkin(8) - t125 * t118;
	t1 = [(-t118 * t130 + t131) * t121 + t122 * t126, (t118 * t129 - t126) * t119 + t121 * t132, t118 * t128 + t127, t114 * t121 + t116 * t124 + 0; -t115 * t122 - t124 * t131, t115 * t119 - t124 * t132, -t118 * t126 + t129, -t114 * t124 + t116 * t121 + 0; t117 * t130 + t118 * t119, t118 * t122 - t120 * t131, -t117 * t123, t118 * pkin(8) + t125 * t117 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:08
	% EndTime: 2020-11-04 22:48:08
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (57->34), mult. (112->60), div. (0->0), fcn. (146->10), ass. (0->29)
	t140 = sin(pkin(6));
	t147 = cos(qJ(3));
	t161 = t140 * t147;
	t142 = sin(qJ(4));
	t148 = cos(qJ(2));
	t160 = t142 * t148;
	t143 = sin(qJ(3));
	t159 = t143 * t140;
	t144 = sin(qJ(2));
	t158 = t144 * t147;
	t145 = sin(qJ(1));
	t157 = t145 * t148;
	t146 = cos(qJ(4));
	t156 = t146 * t148;
	t155 = t147 * t148;
	t149 = cos(qJ(1));
	t154 = t149 * t144;
	t153 = t149 * t148;
	t138 = t147 * pkin(3) + t143 * pkin(10) + pkin(2);
	t152 = t148 * pkin(9) - t138 * t144;
	t151 = t143 * pkin(3) - t147 * pkin(10) + pkin(8);
	t141 = cos(pkin(6));
	t136 = t141 * t158 - t159;
	t150 = t136 * t142 + t141 * t156;
	t137 = t142 * t155 - t146 * t144;
	t135 = t140 * t158 + t141 * t143;
	t134 = t144 * pkin(9) + t138 * t148 + pkin(1);
	t133 = t140 * t151 + t152 * t141;
	t1 = [(-t136 * t145 + t147 * t153) * t146 + (t141 * t157 + t154) * t142, -t149 * t137 + t150 * t145, (-t145 * t141 * t144 + t153) * t143 - t145 * t161, t133 * t145 + t134 * t149 + 0; (t136 * t146 - t141 * t160) * t149 + t145 * (t144 * t142 + t146 * t155), -t145 * t137 - t150 * t149, (t141 * t154 + t157) * t143 + t149 * t161, -t133 * t149 + t134 * t145 + 0; t135 * t146 - t140 * t160, -t135 * t142 - t140 * t156, -t141 * t147 + t144 * t159, -t152 * t140 + t151 * t141 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:08
	% EndTime: 2020-11-04 22:48:08
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (90->36), mult. (125->55), div. (0->0), fcn. (159->12), ass. (0->34)
	t175 = sin(pkin(6));
	t180 = cos(qJ(3));
	t195 = t175 * t180;
	t181 = cos(qJ(2));
	t194 = t175 * t181;
	t177 = sin(qJ(3));
	t193 = t177 * t175;
	t178 = sin(qJ(2));
	t192 = t178 * t180;
	t179 = sin(qJ(1));
	t191 = t179 * t178;
	t190 = t179 * t181;
	t182 = cos(qJ(1));
	t189 = t182 * t178;
	t188 = t182 * t181;
	t171 = cos(qJ(4)) * pkin(4) + pkin(3);
	t183 = pkin(11) + pkin(10);
	t164 = t171 * t180 + t183 * t177 + pkin(2);
	t170 = sin(qJ(4)) * pkin(4) + pkin(9);
	t187 = t164 * t178 - t170 * t181;
	t186 = t171 * t177 - t183 * t180 + pkin(8);
	t176 = cos(pkin(6));
	t166 = t176 * t192 - t193;
	t185 = -t166 * t179 + t180 * t188;
	t184 = t166 * t182 + t180 * t190;
	t174 = qJ(4) + qJ(5);
	t173 = cos(t174);
	t172 = sin(t174);
	t168 = t176 * t190 + t189;
	t167 = t176 * t188 - t191;
	t165 = t175 * t192 + t176 * t177;
	t163 = t164 * t181 + t170 * t178 + pkin(1);
	t162 = t175 * t186 - t187 * t176;
	t1 = [t168 * t172 + t185 * t173, t168 * t173 - t185 * t172, (-t176 * t191 + t188) * t177 - t179 * t195, t162 * t179 + t163 * t182 + 0; -t172 * t167 + t184 * t173, -t173 * t167 - t184 * t172, (t176 * t189 + t190) * t177 + t182 * t195, -t162 * t182 + t163 * t179 + 0; t165 * t173 - t172 * t194, -t165 * t172 - t173 * t194, -t176 * t180 + t178 * t193, t187 * t175 + t186 * t176 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:08
	% EndTime: 2020-11-04 22:48:08
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (107->46), mult. (152->62), div. (0->0), fcn. (198->14), ass. (0->36)
	t217 = sin(qJ(2));
	t230 = pkin(2) * t217;
	t214 = sin(pkin(6));
	t219 = cos(qJ(3));
	t229 = t214 * t219;
	t220 = cos(qJ(2));
	t228 = t214 * t220;
	t216 = sin(qJ(3));
	t227 = t216 * t214;
	t226 = t217 * t219;
	t218 = sin(qJ(1));
	t225 = t218 * t217;
	t224 = t218 * t220;
	t221 = cos(qJ(1));
	t223 = t221 * t217;
	t222 = t221 * t220;
	t213 = qJ(4) + qJ(5);
	t215 = cos(pkin(6));
	t212 = -pkin(12) - pkin(11) - pkin(10);
	t211 = qJ(6) + t213;
	t210 = cos(t211);
	t209 = sin(t211);
	t208 = t220 * pkin(2) + t217 * pkin(9) + pkin(1);
	t207 = pkin(5) * sin(t213) + sin(qJ(4)) * pkin(4);
	t206 = pkin(5) * cos(t213) + cos(qJ(4)) * pkin(4) + pkin(3);
	t205 = t215 * t224 + t223;
	t204 = -t215 * t222 + t225;
	t203 = -t215 * t223 - t224;
	t202 = t214 * t226 + t215 * t216;
	t201 = -t215 * t219 + t217 * t227;
	t200 = t214 * pkin(8) + (pkin(9) * t220 - t230) * t215;
	t199 = -t203 * t219 - t221 * t227;
	t198 = t203 * t216 - t221 * t229;
	t197 = (-t215 * t226 + t227) * t218 + t219 * t222;
	t196 = (t215 * t225 - t222) * t216 + t218 * t229;
	t1 = [t197 * t210 + t205 * t209, -t197 * t209 + t205 * t210, -t196, t196 * t212 + t197 * t206 + t200 * t218 + t205 * t207 + t208 * t221 + 0; t199 * t210 + t204 * t209, -t199 * t209 + t204 * t210, -t198, t198 * t212 + t199 * t206 - t200 * t221 + t204 * t207 + t208 * t218 + 0; t202 * t210 - t209 * t228, -t202 * t209 - t210 * t228, t201, t215 * pkin(8) - t201 * t212 + t202 * t206 + pkin(7) + 0 + (t230 + (-pkin(9) - t207) * t220) * t214; 0, 0, 0, 1;];
	Tc_mdh = t1;
end