% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:46
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:42
	% EndTime: 2019-10-09 22:46:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:42
	% EndTime: 2019-10-09 22:46:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:42
	% EndTime: 2019-10-09 22:46:42
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(11));
	t20 = sin(pkin(6));
	t19 = sin(pkin(11));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:43
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (19->13), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->19)
	t66 = sin(pkin(6));
	t69 = sin(qJ(3));
	t78 = t66 * t69;
	t70 = sin(qJ(2));
	t77 = t66 * t70;
	t71 = cos(qJ(3));
	t76 = t66 * t71;
	t72 = cos(qJ(2));
	t75 = t66 * t72;
	t68 = cos(pkin(6));
	t74 = t68 * t70;
	t73 = t68 * t72;
	t67 = cos(pkin(11));
	t65 = sin(pkin(11));
	t64 = -t65 * t74 + t67 * t72;
	t63 = -t65 * t73 - t67 * t70;
	t62 = t65 * t72 + t67 * t74;
	t61 = -t65 * t70 + t67 * t73;
	t1 = [0, t63 * t71, -t64 * t69 + t65 * t76, 0, 0, 0; 0, t61 * t71, -t62 * t69 - t67 * t76, 0, 0, 0; 0, t71 * t75, t68 * t71 - t69 * t77, 0, 0, 0; 0, -t63 * t69, -t64 * t71 - t65 * t78, 0, 0, 0; 0, -t61 * t69, -t62 * t71 + t67 * t78, 0, 0, 0; 0, -t69 * t75, -t68 * t69 - t70 * t76, 0, 0, 0; 0, t64, 0, 0, 0, 0; 0, t62, 0, 0, 0, 0; 0, t77, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:43
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (59->14), mult. (87->32), div. (0->0), fcn. (134->8), ass. (0->26)
	t93 = sin(pkin(11));
	t94 = sin(pkin(6));
	t104 = t93 * t94;
	t95 = cos(pkin(11));
	t103 = t94 * t95;
	t97 = sin(qJ(2));
	t102 = t94 * t97;
	t98 = cos(qJ(2));
	t101 = t94 * t98;
	t96 = cos(pkin(6));
	t100 = t96 * t97;
	t99 = t96 * t98;
	t92 = qJ(3) + qJ(4);
	t91 = cos(t92);
	t90 = sin(t92);
	t89 = -t93 * t100 + t95 * t98;
	t88 = -t93 * t99 - t95 * t97;
	t87 = t95 * t100 + t93 * t98;
	t86 = -t93 * t97 + t95 * t99;
	t85 = -t91 * t102 - t96 * t90;
	t84 = -t90 * t102 + t96 * t91;
	t83 = -t90 * t104 - t89 * t91;
	t82 = t91 * t104 - t89 * t90;
	t81 = t90 * t103 - t87 * t91;
	t80 = -t91 * t103 - t87 * t90;
	t1 = [0, t88 * t91, t82, t82, 0, 0; 0, t86 * t91, t80, t80, 0, 0; 0, t91 * t101, t84, t84, 0, 0; 0, -t88 * t90, t83, t83, 0, 0; 0, -t86 * t90, t81, t81, 0, 0; 0, -t90 * t101, t85, t85, 0, 0; 0, t89, 0, 0, 0, 0; 0, t87, 0, 0, 0, 0; 0, t102, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:43
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (89->14), mult. (87->32), div. (0->0), fcn. (134->8), ass. (0->26)
	t98 = sin(pkin(11));
	t99 = sin(pkin(6));
	t109 = t98 * t99;
	t100 = cos(pkin(11));
	t108 = t100 * t99;
	t103 = cos(qJ(2));
	t107 = t103 * t99;
	t102 = sin(qJ(2));
	t106 = t99 * t102;
	t101 = cos(pkin(6));
	t105 = t101 * t102;
	t104 = t101 * t103;
	t97 = qJ(3) + qJ(4) + pkin(12);
	t96 = cos(t97);
	t95 = sin(t97);
	t94 = t100 * t103 - t98 * t105;
	t93 = -t100 * t102 - t98 * t104;
	t92 = t100 * t105 + t98 * t103;
	t91 = t100 * t104 - t98 * t102;
	t90 = -t101 * t95 - t96 * t106;
	t89 = t101 * t96 - t95 * t106;
	t88 = -t95 * t109 - t94 * t96;
	t87 = t96 * t109 - t94 * t95;
	t86 = t95 * t108 - t92 * t96;
	t85 = -t96 * t108 - t92 * t95;
	t1 = [0, t93 * t96, t87, t87, 0, 0; 0, t91 * t96, t85, t85, 0, 0; 0, t96 * t107, t89, t89, 0, 0; 0, -t93 * t95, t88, t88, 0, 0; 0, -t91 * t95, t86, t86, 0, 0; 0, -t95 * t107, t90, t90, 0, 0; 0, t94, 0, 0, 0, 0; 0, t92, 0, 0, 0, 0; 0, t106, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:43
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (186->31), mult. (214->65), div. (0->0), fcn. (313->10), ass. (0->37)
	t149 = sin(pkin(11));
	t151 = cos(pkin(11));
	t156 = cos(qJ(2));
	t152 = cos(pkin(6));
	t154 = sin(qJ(2));
	t160 = t152 * t154;
	t142 = t149 * t156 + t151 * t160;
	t148 = qJ(3) + qJ(4) + pkin(12);
	t146 = sin(t148);
	t147 = cos(t148);
	t150 = sin(pkin(6));
	t162 = t150 * t151;
	t134 = -t142 * t146 - t147 * t162;
	t153 = sin(qJ(6));
	t168 = t134 * t153;
	t144 = -t149 * t160 + t151 * t156;
	t163 = t149 * t150;
	t136 = -t144 * t146 + t147 * t163;
	t167 = t136 * t153;
	t161 = t150 * t154;
	t139 = -t146 * t161 + t152 * t147;
	t166 = t139 * t153;
	t165 = t147 * t153;
	t155 = cos(qJ(6));
	t164 = t147 * t155;
	t159 = t152 * t156;
	t158 = t153 * t156;
	t157 = t155 * t156;
	t143 = t149 * t159 + t151 * t154;
	t141 = t149 * t154 - t151 * t159;
	t140 = t152 * t146 + t147 * t161;
	t138 = t139 * t155;
	t137 = t144 * t147 + t146 * t163;
	t135 = t142 * t147 - t146 * t162;
	t133 = t136 * t155;
	t132 = t134 * t155;
	t1 = [0, -t143 * t164 + t144 * t153, t133, t133, 0, -t137 * t153 + t143 * t155; 0, -t141 * t164 + t142 * t153, t132, t132, 0, -t135 * t153 + t141 * t155; 0, (t147 * t157 + t153 * t154) * t150, t138, t138, 0, -t140 * t153 - t150 * t157; 0, t143 * t165 + t144 * t155, -t167, -t167, 0, -t137 * t155 - t143 * t153; 0, t141 * t165 + t142 * t155, -t168, -t168, 0, -t135 * t155 - t141 * t153; 0, (-t147 * t158 + t154 * t155) * t150, -t166, -t166, 0, -t140 * t155 + t150 * t158; 0, -t143 * t146, t137, t137, 0, 0; 0, -t141 * t146, t135, t135, 0, 0; 0, t150 * t156 * t146, t140, t140, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end