% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: S5PRRRR11_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% Jg_rot [3x5]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S5PRRRR11_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_jacobig_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR11_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_jacobig_rot_sym_varpar: pkin has to be [10x1] (double)');
Jg_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	% Symbolic code from jacobig_rot_1_floatb_twist_matlab.m not found
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 1, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (3->2), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->3)
	t44 = sin(pkin(5));
	t43 = pkin(10) + qJ(2);
	t1 = [0, 0, sin(t43) * t44, 0, 0; 0, 0, -cos(t43) * t44, 0, 0; 0, 1, cos(pkin(5)), 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (6->3), mult. (4->2), div. (0->0), fcn. (10->4), ass. (0->6)
	t72 = pkin(10) + qJ(2);
	t73 = sin(pkin(5));
	t75 = cos(t72) * t73;
	t74 = cos(pkin(5));
	t70 = sin(t72) * t73;
	t1 = [0, 0, t70, t70, 0; 0, 0, -t75, -t75, 0; 0, 1, t74, t74, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:37
	% EndTime: 2024-09-27 21:46:37
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (9->4), mult. (6->2), div. (0->0), fcn. (15->4), ass. (0->6)
	t82 = pkin(10) + qJ(2);
	t83 = sin(pkin(5));
	t85 = cos(t82) * t83;
	t84 = cos(pkin(5));
	t80 = sin(t82) * t83;
	t1 = [0, 0, t80, t80, t80; 0, 0, -t85, -t85, -t85; 0, 1, t84, t84, t84;];
	Jg_rot = t1;
end