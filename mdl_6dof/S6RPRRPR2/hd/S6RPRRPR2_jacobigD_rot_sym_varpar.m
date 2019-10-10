% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR2
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRPR2_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR2_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->1), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->2)
	t20 = qJ(1) + pkin(10);
	t1 = [0, 0, qJD(1) * cos(t20), 0, 0, 0; 0, 0, qJD(1) * sin(t20), 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->3), mult. (11->9), div. (0->0), fcn. (11->4), ass. (0->7)
	t77 = sin(qJ(3));
	t80 = qJD(1) * t77;
	t79 = qJD(3) * cos(qJ(3));
	t76 = qJ(1) + pkin(10);
	t75 = cos(t76);
	t74 = sin(t76);
	t1 = [0, 0, qJD(1) * t75, -t74 * t80 + t75 * t79, 0, 0; 0, 0, qJD(1) * t74, t74 * t79 + t75 * t80, 0, 0; 0, 0, 0, qJD(3) * t77, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->3), mult. (11->9), div. (0->0), fcn. (11->4), ass. (0->7)
	t86 = sin(qJ(3));
	t89 = qJD(1) * t86;
	t88 = qJD(3) * cos(qJ(3));
	t85 = qJ(1) + pkin(10);
	t84 = cos(t85);
	t83 = sin(t85);
	t1 = [0, 0, qJD(1) * t84, -t83 * t89 + t84 * t88, 0, 0; 0, 0, qJD(1) * t83, t83 * t88 + t84 * t89, 0, 0; 0, 0, 0, qJD(3) * t86, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (14->3), mult. (20->9), div. (0->0), fcn. (20->4), ass. (0->10)
	t101 = sin(qJ(3));
	t104 = qJD(1) * t101;
	t103 = qJD(3) * cos(qJ(3));
	t100 = qJ(1) + pkin(10);
	t99 = qJD(3) * t101;
	t98 = cos(t100);
	t97 = sin(t100);
	t96 = t97 * t103 + t98 * t104;
	t95 = t98 * t103 - t97 * t104;
	t1 = [0, 0, qJD(1) * t98, t95, 0, t95; 0, 0, qJD(1) * t97, t96, 0, t96; 0, 0, 0, t99, 0, t99;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end