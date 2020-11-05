% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRR12 (for one body)
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
% Datum: 2020-11-04 22:33
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRR12_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:22
	% EndTime: 2020-11-04 22:33:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:22
	% EndTime: 2020-11-04 22:33:22
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t102 = cos(qJ(1));
	t101 = sin(qJ(1));
	t1 = [t102, -t101, 0, 0; t101, t102, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:22
	% EndTime: 2020-11-04 22:33:22
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t103 = sin(pkin(6));
	t106 = sin(qJ(1));
	t114 = t106 * t103;
	t105 = sin(qJ(2));
	t113 = t106 * t105;
	t107 = cos(qJ(2));
	t112 = t106 * t107;
	t108 = cos(qJ(1));
	t111 = t108 * t103;
	t110 = t108 * t105;
	t109 = t108 * t107;
	t104 = cos(pkin(6));
	t1 = [-t104 * t113 + t109, -t104 * t112 - t110, t114, t108 * pkin(1) + pkin(8) * t114 + 0; t104 * t110 + t112, t104 * t109 - t113, -t111, t106 * pkin(1) - pkin(8) * t111 + 0; t103 * t105, t103 * t107, t104, t104 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:22
	% EndTime: 2020-11-04 22:33:22
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t118 = sin(pkin(6));
	t123 = cos(qJ(3));
	t133 = t118 * t123;
	t120 = sin(qJ(3));
	t132 = t120 * t118;
	t121 = sin(qJ(2));
	t131 = t121 * t123;
	t122 = sin(qJ(1));
	t130 = t122 * t121;
	t124 = cos(qJ(2));
	t129 = t122 * t124;
	t125 = cos(qJ(1));
	t128 = t125 * t121;
	t127 = t125 * t124;
	t126 = pkin(2) * t121 - pkin(9) * t124;
	t119 = cos(pkin(6));
	t117 = t124 * pkin(2) + t121 * pkin(9) + pkin(1);
	t116 = t119 * t128 + t129;
	t115 = t118 * pkin(8) - t126 * t119;
	t1 = [(-t119 * t131 + t132) * t122 + t123 * t127, (t119 * t130 - t127) * t120 + t122 * t133, t119 * t129 + t128, t115 * t122 + t117 * t125 + 0; t116 * t123 - t125 * t132, -t116 * t120 - t125 * t133, -t119 * t127 + t130, -t115 * t125 + t117 * t122 + 0; t118 * t131 + t119 * t120, t119 * t123 - t121 * t132, -t118 * t124, t119 * pkin(8) + t126 * t118 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:22
	% EndTime: 2020-11-04 22:33:22
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (57->32), mult. (112->56), div. (0->0), fcn. (146->10), ass. (0->29)
	t144 = sin(pkin(6));
	t150 = cos(qJ(3));
	t161 = t144 * t150;
	t151 = cos(qJ(2));
	t160 = t144 * t151;
	t146 = cos(pkin(6));
	t148 = sin(qJ(2));
	t159 = t146 * t148;
	t158 = t146 * t151;
	t147 = sin(qJ(3));
	t157 = t147 * t144;
	t156 = t148 * t150;
	t155 = t150 * t151;
	t142 = pkin(3) * t150 + qJ(4) * t147 + pkin(2);
	t154 = t151 * pkin(9) - t142 * t148;
	t153 = t147 * pkin(3) - qJ(4) * t150 + pkin(8);
	t152 = cos(qJ(1));
	t149 = sin(qJ(1));
	t145 = cos(pkin(12));
	t143 = sin(pkin(12));
	t141 = t143 * t148 + t145 * t155;
	t140 = t146 * t156 - t157;
	t139 = t144 * t156 + t146 * t147;
	t138 = t143 * t155 - t145 * t148;
	t137 = t148 * pkin(9) + t142 * t151 + pkin(1);
	t136 = t143 * t140 + t145 * t158;
	t135 = -t145 * t140 + t143 * t158;
	t134 = t144 * t153 + t154 * t146;
	t1 = [t135 * t149 + t152 * t141, t136 * t149 - t152 * t138, (-t149 * t159 + t152 * t151) * t147 - t149 * t161, t134 * t149 + t137 * t152 + 0; -t135 * t152 + t149 * t141, -t136 * t152 - t149 * t138, (t149 * t151 + t152 * t159) * t147 + t152 * t161, -t134 * t152 + t137 * t149 + 0; t139 * t145 - t143 * t160, -t139 * t143 - t145 * t160, -t146 * t150 + t148 * t157, -t154 * t144 + t153 * t146 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:22
	% EndTime: 2020-11-04 22:33:22
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (90->36), mult. (125->55), div. (0->0), fcn. (159->12), ass. (0->34)
	t175 = sin(pkin(6));
	t181 = cos(qJ(3));
	t195 = t175 * t181;
	t182 = cos(qJ(2));
	t194 = t175 * t182;
	t178 = sin(qJ(3));
	t193 = t178 * t175;
	t179 = sin(qJ(2));
	t192 = t179 * t181;
	t183 = cos(qJ(1));
	t191 = t179 * t183;
	t180 = sin(qJ(1));
	t190 = t180 * t179;
	t189 = t180 * t182;
	t188 = t182 * t183;
	t171 = cos(pkin(12)) * pkin(4) + pkin(3);
	t177 = qJ(4) + pkin(10);
	t164 = t171 * t181 + t177 * t178 + pkin(2);
	t170 = sin(pkin(12)) * pkin(4) + pkin(9);
	t187 = t164 * t179 - t170 * t182;
	t186 = t171 * t178 - t177 * t181 + pkin(8);
	t176 = cos(pkin(6));
	t166 = t176 * t192 - t193;
	t185 = -t166 * t180 + t181 * t188;
	t184 = t166 * t183 + t181 * t189;
	t174 = pkin(12) + qJ(5);
	t173 = cos(t174);
	t172 = sin(t174);
	t168 = t176 * t189 + t191;
	t167 = t176 * t188 - t190;
	t165 = t175 * t192 + t176 * t178;
	t163 = t164 * t182 + t170 * t179 + pkin(1);
	t162 = t186 * t175 - t187 * t176;
	t1 = [t168 * t172 + t185 * t173, t168 * t173 - t185 * t172, (-t176 * t190 + t188) * t178 - t180 * t195, t162 * t180 + t163 * t183 + 0; -t172 * t167 + t184 * t173, -t173 * t167 - t184 * t172, (t176 * t191 + t189) * t178 + t183 * t195, -t162 * t183 + t163 * t180 + 0; t165 * t173 - t172 * t194, -t165 * t172 - t173 * t194, -t176 * t181 + t179 * t193, t187 * t175 + t186 * t176 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:22
	% EndTime: 2020-11-04 22:33:22
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
	t213 = pkin(12) + qJ(5);
	t215 = cos(pkin(6));
	t212 = -pkin(11) - pkin(10) - qJ(4);
	t211 = qJ(6) + t213;
	t210 = cos(t211);
	t209 = sin(t211);
	t208 = t220 * pkin(2) + t217 * pkin(9) + pkin(1);
	t207 = pkin(5) * sin(t213) + sin(pkin(12)) * pkin(4);
	t206 = pkin(5) * cos(t213) + cos(pkin(12)) * pkin(4) + pkin(3);
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