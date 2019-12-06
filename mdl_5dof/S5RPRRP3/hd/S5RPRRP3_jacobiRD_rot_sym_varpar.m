% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRRP3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:04:27
	% EndTime: 2019-12-05 18:04:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:04:28
	% EndTime: 2019-12-05 18:04:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = qJD(1) * cos(qJ(1));
	t7 = qJD(1) * sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t7, 0, 0, 0, 0; -t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; t9, 0, 0, 0, 0; t7, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:04:28
	% EndTime: 2019-12-05 18:04:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->2), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(8);
	t13 = qJD(1) * cos(t12);
	t10 = qJD(1) * sin(t12);
	t1 = [0, 0, 0, 0, 0; t10, 0, 0, 0, 0; -t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; t13, 0, 0, 0, 0; t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:04:28
	% EndTime: 2019-12-05 18:04:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->11), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t85 = sin(qJ(3));
	t90 = qJD(1) * t85;
	t86 = cos(qJ(3));
	t89 = qJD(1) * t86;
	t88 = qJD(3) * t85;
	t87 = qJD(3) * t86;
	t84 = qJ(1) + pkin(8);
	t83 = cos(t84);
	t82 = sin(t84);
	t81 = -t82 * t88 + t83 * t89;
	t80 = t82 * t87 + t83 * t90;
	t79 = t82 * t89 + t83 * t88;
	t78 = t82 * t90 - t83 * t87;
	t1 = [0, 0, -t88, 0, 0; t79, 0, t80, 0, 0; -t81, 0, t78, 0, 0; 0, 0, -t87, 0, 0; -t78, 0, t81, 0, 0; t80, 0, t79, 0, 0; 0, 0, 0, 0, 0; -qJD(1) * t83, 0, 0, 0, 0; -qJD(1) * t82, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:04:28
	% EndTime: 2019-12-05 18:04:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (88->15), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->16)
	t140 = qJ(3) + qJ(4);
	t136 = sin(t140);
	t138 = qJD(3) + qJD(4);
	t144 = t138 * t136;
	t137 = cos(t140);
	t143 = t138 * t137;
	t142 = qJD(1) * t136;
	t141 = qJD(1) * t137;
	t139 = qJ(1) + pkin(8);
	t135 = cos(t139);
	t134 = sin(t139);
	t133 = -t134 * t144 + t135 * t141;
	t132 = t134 * t143 + t135 * t142;
	t131 = t134 * t141 + t135 * t144;
	t130 = t134 * t142 - t135 * t143;
	t1 = [0, 0, -t144, -t144, 0; t131, 0, t132, t132, 0; -t133, 0, t130, t130, 0; 0, 0, -t143, -t143, 0; -t130, 0, t133, t133, 0; t132, 0, t131, t131, 0; 0, 0, 0, 0, 0; -qJD(1) * t135, 0, 0, 0, 0; -qJD(1) * t134, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:04:28
	% EndTime: 2019-12-05 18:04:28
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (88->15), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->16)
	t142 = qJ(3) + qJ(4);
	t138 = sin(t142);
	t140 = qJD(3) + qJD(4);
	t146 = t140 * t138;
	t139 = cos(t142);
	t145 = t140 * t139;
	t144 = qJD(1) * t138;
	t143 = qJD(1) * t139;
	t141 = qJ(1) + pkin(8);
	t137 = cos(t141);
	t136 = sin(t141);
	t135 = -t136 * t146 + t137 * t143;
	t134 = t136 * t145 + t137 * t144;
	t133 = t136 * t143 + t137 * t146;
	t132 = t136 * t144 - t137 * t145;
	t1 = [0, 0, -t146, -t146, 0; t133, 0, t134, t134, 0; -t135, 0, t132, t132, 0; 0, 0, -t145, -t145, 0; -t132, 0, t135, t135, 0; t134, 0, t133, t133, 0; 0, 0, 0, 0, 0; -qJD(1) * t137, 0, 0, 0, 0; -qJD(1) * t136, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end