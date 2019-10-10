% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:24
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPPRR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(10));
	t20 = sin(pkin(6));
	t19 = sin(pkin(10));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (14->7), mult. (40->16), div. (0->0), fcn. (60->8), ass. (0->13)
	t35 = sin(pkin(11));
	t38 = cos(pkin(11));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t34 = -t42 * t35 - t41 * t38;
	t33 = t41 * t35 - t42 * t38;
	t40 = cos(pkin(6));
	t39 = cos(pkin(10));
	t37 = sin(pkin(6));
	t36 = sin(pkin(10));
	t32 = t34 * t40;
	t31 = t33 * t40;
	t1 = [0, t36 * t31 + t39 * t34, 0, 0, 0, 0; 0, -t39 * t31 + t36 * t34, 0, 0, 0, 0; 0, -t33 * t37, 0, 0, 0, 0; 0, -t36 * t32 + t39 * t33, 0, 0, 0, 0; 0, t39 * t32 + t36 * t33, 0, 0, 0, 0; 0, t34 * t37, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:28
	% EndTime: 2019-10-09 21:24:28
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (24->9), mult. (66->22), div. (0->0), fcn. (96->10), ass. (0->18)
	t74 = sin(pkin(11));
	t78 = cos(pkin(11));
	t81 = sin(qJ(2));
	t82 = cos(qJ(2));
	t83 = t82 * t74 + t81 * t78;
	t71 = t81 * t74 - t82 * t78;
	t80 = cos(pkin(6));
	t79 = cos(pkin(10));
	t77 = cos(pkin(12));
	t76 = sin(pkin(6));
	t75 = sin(pkin(10));
	t73 = sin(pkin(12));
	t70 = t83 * t80;
	t69 = t71 * t80;
	t68 = t71 * t76;
	t67 = t75 * t69 - t79 * t83;
	t66 = -t79 * t69 - t75 * t83;
	t1 = [0, t67 * t77, 0, 0, 0, 0; 0, t66 * t77, 0, 0, 0, 0; 0, -t68 * t77, 0, 0, 0, 0; 0, -t67 * t73, 0, 0, 0, 0; 0, -t66 * t73, 0, 0, 0, 0; 0, t68 * t73, 0, 0, 0, 0; 0, -t75 * t70 - t79 * t71, 0, 0, 0, 0; 0, t79 * t70 - t75 * t71, 0, 0, 0, 0; 0, t83 * t76, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:28
	% EndTime: 2019-10-09 21:24:28
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (62->16), mult. (122->36), div. (0->0), fcn. (178->10), ass. (0->24)
	t102 = sin(pkin(10));
	t103 = sin(pkin(6));
	t111 = t102 * t103;
	t105 = cos(pkin(10));
	t110 = t103 * t105;
	t106 = cos(pkin(6));
	t101 = sin(pkin(11));
	t104 = cos(pkin(11));
	t107 = sin(qJ(2));
	t108 = cos(qJ(2));
	t109 = t108 * t101 + t107 * t104;
	t94 = t109 * t106;
	t95 = t107 * t101 - t108 * t104;
	t88 = -t102 * t95 + t105 * t94;
	t90 = -t102 * t94 - t105 * t95;
	t100 = pkin(12) + qJ(5);
	t99 = cos(t100);
	t98 = sin(t100);
	t93 = t95 * t106;
	t92 = t109 * t103;
	t91 = t95 * t103;
	t89 = t102 * t93 - t105 * t109;
	t87 = -t102 * t109 - t105 * t93;
	t1 = [0, t89 * t99, 0, 0, t111 * t99 - t90 * t98, 0; 0, t87 * t99, 0, 0, -t110 * t99 - t88 * t98, 0; 0, -t91 * t99, 0, 0, t106 * t99 - t92 * t98, 0; 0, -t89 * t98, 0, 0, -t111 * t98 - t90 * t99, 0; 0, -t87 * t98, 0, 0, t110 * t98 - t88 * t99, 0; 0, t91 * t98, 0, 0, -t106 * t98 - t92 * t99, 0; 0, t90, 0, 0, 0, 0; 0, t88, 0, 0, 0, 0; 0, t92, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:28
	% EndTime: 2019-10-09 21:24:28
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (153->29), mult. (319->65), div. (0->0), fcn. (454->12), ass. (0->34)
	t153 = sin(pkin(11));
	t156 = cos(pkin(11));
	t160 = sin(qJ(2));
	t162 = cos(qJ(2));
	t146 = t160 * t153 - t162 * t156;
	t152 = pkin(12) + qJ(5);
	t151 = cos(t152);
	t159 = sin(qJ(6));
	t172 = t151 * t159;
	t161 = cos(qJ(6));
	t171 = t151 * t161;
	t154 = sin(pkin(10));
	t155 = sin(pkin(6));
	t170 = t154 * t155;
	t157 = cos(pkin(10));
	t169 = t155 * t157;
	t158 = cos(pkin(6));
	t164 = t162 * t153 + t160 * t156;
	t145 = t164 * t158;
	t166 = t157 * t145 - t154 * t146;
	t165 = -t154 * t145 - t157 * t146;
	t163 = t146 * t158;
	t150 = sin(t152);
	t144 = t164 * t155;
	t143 = t146 * t155;
	t141 = t144 * t151 + t158 * t150;
	t140 = -t144 * t150 + t158 * t151;
	t138 = t154 * t163 - t157 * t164;
	t135 = -t154 * t164 - t157 * t163;
	t133 = t150 * t170 + t151 * t165;
	t132 = -t150 * t165 + t151 * t170;
	t131 = -t150 * t169 + t151 * t166;
	t130 = -t150 * t166 - t151 * t169;
	t1 = [0, t138 * t171 + t159 * t165, 0, 0, t132 * t161, -t133 * t159 - t138 * t161; 0, t135 * t171 + t159 * t166, 0, 0, t130 * t161, -t131 * t159 - t135 * t161; 0, -t143 * t171 + t144 * t159, 0, 0, t140 * t161, -t141 * t159 + t143 * t161; 0, -t138 * t172 + t161 * t165, 0, 0, -t132 * t159, -t133 * t161 + t138 * t159; 0, -t135 * t172 + t161 * t166, 0, 0, -t130 * t159, -t131 * t161 + t135 * t159; 0, t143 * t172 + t144 * t161, 0, 0, -t140 * t159, -t141 * t161 - t143 * t159; 0, t138 * t150, 0, 0, t133, 0; 0, t135 * t150, 0, 0, t131, 0; 0, -t143 * t150, 0, 0, t141, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end