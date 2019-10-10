% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:59
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRR4_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR4_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:59:27
	% EndTime: 2019-10-09 21:59:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:59:27
	% EndTime: 2019-10-09 21:59:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:59:27
	% EndTime: 2019-10-09 21:59:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(11)) * t18, 0, 0, 0, 0; 0, -cos(pkin(11)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:59:27
	% EndTime: 2019-10-09 21:59:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t50 = sin(pkin(6));
	t1 = [0, sin(pkin(11)) * t50, 0, 0, 0, 0; 0, -cos(pkin(11)) * t50, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:59:27
	% EndTime: 2019-10-09 21:59:27
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t61 = cos(pkin(6));
	t63 = cos(qJ(2));
	t64 = t61 * t63;
	t62 = sin(qJ(2));
	t60 = cos(pkin(11));
	t59 = sin(pkin(6));
	t58 = sin(pkin(11));
	t1 = [0, t58 * t59, 0, t58 * t64 + t60 * t62, 0, 0; 0, -t60 * t59, 0, t58 * t62 - t60 * t64, 0, 0; 0, t61, 0, -t59 * t63, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:59:27
	% EndTime: 2019-10-09 21:59:27
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (15->10), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->14)
	t92 = sin(pkin(11));
	t93 = sin(pkin(6));
	t101 = t92 * t93;
	t94 = cos(pkin(11));
	t100 = t94 * t93;
	t95 = cos(pkin(6));
	t96 = sin(qJ(2));
	t99 = t95 * t96;
	t97 = cos(qJ(2));
	t98 = t95 * t97;
	t91 = pkin(12) + qJ(4);
	t90 = cos(t91);
	t89 = sin(t91);
	t1 = [0, t101, 0, t92 * t98 + t94 * t96, (-t92 * t99 + t94 * t97) * t89 - t90 * t101, 0; 0, -t100, 0, t92 * t96 - t94 * t98, (t92 * t97 + t94 * t99) * t89 + t90 * t100, 0; 0, t95, 0, -t93 * t97, t93 * t96 * t89 - t95 * t90, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:59:27
	% EndTime: 2019-10-09 21:59:27
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (26->10), mult. (39->20), div. (0->0), fcn. (63->8), ass. (0->17)
	t111 = sin(pkin(11));
	t112 = sin(pkin(6));
	t120 = t111 * t112;
	t113 = cos(pkin(11));
	t119 = t113 * t112;
	t114 = cos(pkin(6));
	t115 = sin(qJ(2));
	t118 = t114 * t115;
	t116 = cos(qJ(2));
	t117 = t114 * t116;
	t110 = pkin(12) + qJ(4);
	t109 = cos(t110);
	t108 = sin(t110);
	t107 = t112 * t115 * t108 - t114 * t109;
	t106 = (-t111 * t118 + t113 * t116) * t108 - t109 * t120;
	t105 = (t111 * t116 + t113 * t118) * t108 + t109 * t119;
	t1 = [0, t120, 0, t111 * t117 + t113 * t115, t106, t106; 0, -t119, 0, t111 * t115 - t113 * t117, t105, t105; 0, t114, 0, -t112 * t116, t107, t107;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end