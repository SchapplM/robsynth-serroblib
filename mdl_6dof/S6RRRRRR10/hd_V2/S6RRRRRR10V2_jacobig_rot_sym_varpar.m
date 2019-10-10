% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:38
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR10V2_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR10V2_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobig_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)), 0, 0, 0, 0; 0, -cos(qJ(1)), 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t17 = cos(qJ(1));
	t16 = sin(qJ(1));
	t1 = [0, t16, t16, 0, 0, 0; 0, -t17, -t17, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->4), mult. (2->2), div. (0->0), fcn. (9->4), ass. (0->5)
	t87 = cos(qJ(1));
	t86 = sin(qJ(1));
	t85 = qJ(2) + qJ(3);
	t84 = sin(t85);
	t1 = [0, t86, t86, t87 * t84, 0, 0; 0, -t87, -t87, t86 * t84, 0, 0; 1, 0, 0, -cos(t85), 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->6), mult. (9->8), div. (0->0), fcn. (21->6), ass. (0->9)
	t108 = qJ(2) + qJ(3);
	t107 = cos(t108);
	t109 = sin(qJ(4));
	t113 = t107 * t109;
	t112 = cos(qJ(1));
	t111 = cos(qJ(4));
	t110 = sin(qJ(1));
	t106 = sin(t108);
	t1 = [0, t110, t110, t112 * t106, -t110 * t111 + t112 * t113, 0; 0, -t112, -t112, t110 * t106, t110 * t113 + t112 * t111, 0; 1, 0, 0, -t107, t106 * t109, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (22->11), mult. (24->18), div. (0->0), fcn. (44->8), ass. (0->16)
	t154 = qJ(2) + qJ(3);
	t152 = sin(t154);
	t157 = sin(qJ(1));
	t166 = t157 * t152;
	t156 = sin(qJ(4));
	t165 = t157 * t156;
	t159 = cos(qJ(4));
	t164 = t157 * t159;
	t160 = cos(qJ(1));
	t163 = t160 * t152;
	t162 = t160 * t156;
	t161 = t160 * t159;
	t158 = cos(qJ(5));
	t155 = sin(qJ(5));
	t153 = cos(t154);
	t1 = [0, t157, t157, t163, t153 * t162 - t164, (t153 * t161 + t165) * t155 - t158 * t163; 0, -t160, -t160, t166, t153 * t165 + t161, (t153 * t164 - t162) * t155 - t158 * t166; 1, 0, 0, -t153, t152 * t156, t152 * t159 * t155 + t153 * t158;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end