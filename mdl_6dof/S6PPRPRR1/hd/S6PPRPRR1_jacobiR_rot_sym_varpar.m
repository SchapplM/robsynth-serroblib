% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:08
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRPRR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRPRR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (20->14), mult. (62->33), div. (0->0), fcn. (88->10), ass. (0->20)
	t53 = sin(pkin(11));
	t59 = cos(pkin(6));
	t68 = t53 * t59;
	t54 = sin(pkin(7));
	t55 = sin(pkin(6));
	t67 = t54 * t55;
	t66 = t54 * t59;
	t56 = cos(pkin(12));
	t58 = cos(pkin(7));
	t65 = t56 * t58;
	t57 = cos(pkin(11));
	t64 = t57 * t59;
	t52 = sin(pkin(12));
	t63 = -(-t53 * t52 + t56 * t64) * t58 + t57 * t67;
	t62 = (-t57 * t52 - t56 * t68) * t58 + t53 * t67;
	t61 = cos(qJ(3));
	t60 = sin(qJ(3));
	t51 = -t52 * t68 + t57 * t56;
	t49 = t52 * t64 + t53 * t56;
	t1 = [0, 0, -t51 * t60 + t62 * t61, 0, 0, 0; 0, 0, -t49 * t60 - t63 * t61, 0, 0, 0; 0, 0, t61 * t66 + (-t52 * t60 + t61 * t65) * t55, 0, 0, 0; 0, 0, -t51 * t61 - t62 * t60, 0, 0, 0; 0, 0, -t49 * t61 + t63 * t60, 0, 0, 0; 0, 0, -t60 * t66 + (-t52 * t61 - t60 * t65) * t55, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (38->18), mult. (110->40), div. (0->0), fcn. (154->12), ass. (0->27)
	t75 = sin(pkin(11));
	t77 = sin(pkin(6));
	t88 = t75 * t77;
	t82 = cos(pkin(6));
	t87 = t75 * t82;
	t80 = cos(pkin(11));
	t86 = t77 * t80;
	t85 = t80 * t82;
	t73 = sin(pkin(13));
	t78 = cos(pkin(13));
	t83 = sin(qJ(3));
	t84 = cos(qJ(3));
	t72 = -t84 * t73 - t83 * t78;
	t71 = t83 * t73 - t84 * t78;
	t81 = cos(pkin(7));
	t79 = cos(pkin(12));
	t76 = sin(pkin(7));
	t74 = sin(pkin(12));
	t70 = -t74 * t87 + t80 * t79;
	t69 = -t80 * t74 - t79 * t87;
	t68 = t74 * t85 + t75 * t79;
	t67 = -t75 * t74 + t79 * t85;
	t66 = t72 * t81;
	t65 = t71 * t81;
	t64 = t72 * t76;
	t63 = t71 * t76;
	t1 = [0, 0, -t63 * t88 - t69 * t65 + t70 * t72, 0, 0, 0; 0, 0, t63 * t86 - t67 * t65 + t68 * t72, 0, 0, 0; 0, 0, -t82 * t63 + (-t65 * t79 + t72 * t74) * t77, 0, 0, 0; 0, 0, t64 * t88 + t69 * t66 + t70 * t71, 0, 0, 0; 0, 0, -t64 * t86 + t67 * t66 + t68 * t71, 0, 0, 0; 0, 0, t82 * t64 + (t66 * t79 + t71 * t74) * t77, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:44
	% EndTime: 2019-10-09 21:08:44
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (114->30), mult. (323->66), div. (0->0), fcn. (449->14), ass. (0->39)
	t145 = sin(pkin(11));
	t147 = sin(pkin(6));
	t162 = t145 * t147;
	t152 = cos(pkin(6));
	t161 = t145 * t152;
	t150 = cos(pkin(11));
	t160 = t147 * t150;
	t151 = cos(pkin(7));
	t159 = t147 * t151;
	t158 = t150 * t152;
	t143 = sin(pkin(13));
	t148 = cos(pkin(13));
	t154 = sin(qJ(3));
	t156 = cos(qJ(3));
	t157 = t156 * t143 + t154 * t148;
	t140 = t154 * t143 - t156 * t148;
	t146 = sin(pkin(7));
	t132 = t157 * t146;
	t134 = t157 * t151;
	t144 = sin(pkin(12));
	t149 = cos(pkin(12));
	t136 = -t145 * t144 + t149 * t158;
	t137 = t144 * t158 + t145 * t149;
	t124 = -t132 * t160 + t136 * t134 - t137 * t140;
	t138 = -t150 * t144 - t149 * t161;
	t139 = -t144 * t161 + t150 * t149;
	t126 = t132 * t162 + t138 * t134 - t139 * t140;
	t128 = t152 * t132 + (t134 * t149 - t140 * t144) * t147;
	t155 = cos(qJ(5));
	t153 = sin(qJ(5));
	t135 = -t147 * t149 * t146 + t152 * t151;
	t133 = t140 * t151;
	t131 = t140 * t146;
	t130 = -t138 * t146 + t145 * t159;
	t129 = -t136 * t146 - t150 * t159;
	t127 = -t152 * t131 + (-t133 * t149 - t144 * t157) * t147;
	t125 = -t131 * t162 - t138 * t133 - t139 * t157;
	t123 = t131 * t160 - t136 * t133 - t137 * t157;
	t1 = [0, 0, t125 * t155, 0, -t126 * t153 + t130 * t155, 0; 0, 0, t123 * t155, 0, -t124 * t153 + t129 * t155, 0; 0, 0, t127 * t155, 0, -t128 * t153 + t135 * t155, 0; 0, 0, -t125 * t153, 0, -t126 * t155 - t130 * t153, 0; 0, 0, -t123 * t153, 0, -t124 * t155 - t129 * t153, 0; 0, 0, -t127 * t153, 0, -t128 * t155 - t135 * t153, 0; 0, 0, t126, 0, 0, 0; 0, 0, t124, 0, 0, 0; 0, 0, t128, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:44
	% EndTime: 2019-10-09 21:08:44
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (283->42), mult. (804->95), div. (0->0), fcn. (1108->16), ass. (0->49)
	t203 = sin(pkin(11));
	t205 = sin(pkin(6));
	t227 = t203 * t205;
	t210 = cos(pkin(6));
	t226 = t203 * t210;
	t208 = cos(pkin(11));
	t225 = t205 * t208;
	t209 = cos(pkin(7));
	t224 = t205 * t209;
	t223 = t208 * t210;
	t211 = sin(qJ(6));
	t215 = cos(qJ(5));
	t222 = t211 * t215;
	t214 = cos(qJ(6));
	t221 = t214 * t215;
	t201 = sin(pkin(13));
	t206 = cos(pkin(13));
	t213 = sin(qJ(3));
	t216 = cos(qJ(3));
	t220 = t216 * t201 + t213 * t206;
	t198 = t213 * t201 - t216 * t206;
	t204 = sin(pkin(7));
	t190 = t220 * t204;
	t192 = t220 * t209;
	t202 = sin(pkin(12));
	t207 = cos(pkin(12));
	t194 = -t203 * t202 + t207 * t223;
	t195 = t202 * t223 + t203 * t207;
	t219 = -t190 * t225 + t194 * t192 - t195 * t198;
	t196 = -t208 * t202 - t207 * t226;
	t197 = -t202 * t226 + t208 * t207;
	t218 = t190 * t227 + t196 * t192 - t197 * t198;
	t217 = t210 * t190 + (t192 * t207 - t198 * t202) * t205;
	t212 = sin(qJ(5));
	t193 = -t205 * t207 * t204 + t210 * t209;
	t191 = t198 * t209;
	t189 = t198 * t204;
	t187 = -t196 * t204 + t203 * t224;
	t186 = -t194 * t204 - t208 * t224;
	t184 = -t210 * t189 + (-t191 * t207 - t202 * t220) * t205;
	t182 = t193 * t212 + t215 * t217;
	t181 = t193 * t215 - t212 * t217;
	t179 = -t189 * t227 - t196 * t191 - t197 * t220;
	t176 = t189 * t225 - t194 * t191 - t195 * t220;
	t174 = t187 * t212 + t215 * t218;
	t173 = t187 * t215 - t212 * t218;
	t172 = t186 * t212 + t215 * t219;
	t171 = t186 * t215 - t212 * t219;
	t1 = [0, 0, t179 * t221 + t211 * t218, 0, t173 * t214, -t174 * t211 - t179 * t214; 0, 0, t176 * t221 + t211 * t219, 0, t171 * t214, -t172 * t211 - t176 * t214; 0, 0, t184 * t221 + t211 * t217, 0, t181 * t214, -t182 * t211 - t184 * t214; 0, 0, -t179 * t222 + t214 * t218, 0, -t173 * t211, -t174 * t214 + t179 * t211; 0, 0, -t176 * t222 + t214 * t219, 0, -t171 * t211, -t172 * t214 + t176 * t211; 0, 0, -t184 * t222 + t214 * t217, 0, -t181 * t211, -t182 * t214 + t184 * t211; 0, 0, t179 * t212, 0, t174, 0; 0, 0, t176 * t212, 0, t172, 0; 0, 0, t184 * t212, 0, t182, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end