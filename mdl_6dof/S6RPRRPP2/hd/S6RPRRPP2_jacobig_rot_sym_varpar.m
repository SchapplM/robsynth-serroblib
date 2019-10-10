% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:11
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRPP2_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP2_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_jacobig_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->2)
	t10 = qJ(1) + pkin(9);
	t1 = [0, 0, sin(t10), 0, 0, 0; 0, 0, -cos(t10), 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->3), mult. (2->2), div. (0->0), fcn. (7->4), ass. (0->5)
	t58 = sin(qJ(3));
	t57 = qJ(1) + pkin(9);
	t56 = cos(t57);
	t55 = sin(t57);
	t1 = [0, 0, t55, t56 * t58, 0, 0; 0, 0, -t56, t55 * t58, 0, 0; 1, 0, 0, -cos(qJ(3)), 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->3), mult. (2->2), div. (0->0), fcn. (7->4), ass. (0->5)
	t71 = sin(qJ(3));
	t70 = qJ(1) + pkin(9);
	t69 = cos(t70);
	t68 = sin(t70);
	t1 = [0, 0, t68, t69 * t71, 0, 0; 0, 0, -t69, t68 * t71, 0, 0; 1, 0, 0, -cos(qJ(3)), 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (6->3), mult. (2->2), div. (0->0), fcn. (7->4), ass. (0->5)
	t63 = sin(qJ(3));
	t62 = qJ(1) + pkin(9);
	t61 = cos(t62);
	t60 = sin(t62);
	t1 = [0, 0, t60, t61 * t63, 0, 0; 0, 0, -t61, t60 * t63, 0, 0; 1, 0, 0, -cos(qJ(3)), 0, 0;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end