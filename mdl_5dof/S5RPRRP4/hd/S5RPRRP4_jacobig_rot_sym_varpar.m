% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRP4
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
%   Siehe auch: S5RPRRP4_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% Jg_rot [3x5]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S5RPRRP4_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_jacobig_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP4_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_jacobig_rot_sym_varpar: pkin has to be [8x1] (double)');
Jg_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t44 = sin(pkin(8));
	t1 = [0, 0, cos(qJ(1)) * t44, 0, 0; 0, 0, sin(qJ(1)) * t44, 0, 0; 1, 0, -cos(pkin(8)), 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (4->2), div. (0->0), fcn. (10->4), ass. (0->5)
	t60 = cos(pkin(8));
	t59 = sin(pkin(8));
	t58 = cos(qJ(1)) * t59;
	t57 = sin(qJ(1)) * t59;
	t1 = [0, 0, t58, t58, 0; 0, 0, t57, t57, 0; 1, 0, -t60, -t60, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (4->2), div. (0->0), fcn. (10->4), ass. (0->5)
	t64 = cos(pkin(8));
	t63 = sin(pkin(8));
	t62 = cos(qJ(1)) * t63;
	t61 = sin(qJ(1)) * t63;
	t1 = [0, 0, t62, t62, 0; 0, 0, t61, t61, 0; 1, 0, -t64, -t64, 0;];
	Jg_rot = t1;
end