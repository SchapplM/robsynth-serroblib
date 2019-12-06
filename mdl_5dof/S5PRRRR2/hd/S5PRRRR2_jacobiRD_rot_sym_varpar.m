% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRRR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRRRR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_jacobiRD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:09
	% EndTime: 2019-12-05 17:05:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:09
	% EndTime: 2019-12-05 17:05:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:10
	% EndTime: 2019-12-05 17:05:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(2) * sin(qJ(2));
	t30 = qJD(2) * cos(qJ(2));
	t1 = [0, -t30, 0, 0, 0; 0, -t31, 0, 0, 0; 0, 0, 0, 0, 0; 0, t31, 0, 0, 0; 0, -t30, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:09
	% EndTime: 2019-12-05 17:05:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t47 = qJD(2) + qJD(3);
	t48 = qJ(2) + qJ(3);
	t49 = t47 * cos(t48);
	t44 = t47 * sin(t48);
	t1 = [0, -t49, -t49, 0, 0; 0, -t44, -t44, 0, 0; 0, 0, 0, 0, 0; 0, t44, t44, 0, 0; 0, -t49, -t49, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:09
	% EndTime: 2019-12-05 17:05:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (57->11), mult. (12->2), div. (0->0), fcn. (12->2), ass. (0->5)
	t59 = qJD(2) + qJD(3) + qJD(4);
	t60 = qJ(2) + qJ(3) + qJ(4);
	t61 = t59 * cos(t60);
	t56 = t59 * sin(t60);
	t1 = [0, -t61, -t61, -t61, 0; 0, -t56, -t56, -t56, 0; 0, 0, 0, 0, 0; 0, t56, t56, t56, 0; 0, -t61, -t61, -t61, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:10
	% EndTime: 2019-12-05 17:05:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (141->15), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t80 = qJ(2) + qJ(3) + qJ(4);
	t77 = sin(t80);
	t79 = qJD(2) + qJD(3) + qJD(4);
	t87 = t79 * t77;
	t81 = sin(qJ(5));
	t86 = t79 * t81;
	t82 = cos(qJ(5));
	t85 = t79 * t82;
	t84 = qJD(5) * t81;
	t83 = qJD(5) * t82;
	t78 = cos(t80);
	t76 = t79 * t78;
	t75 = t77 * t84 - t78 * t85;
	t74 = t77 * t83 + t78 * t86;
	t73 = t77 * t85 + t78 * t84;
	t72 = t77 * t86 - t78 * t83;
	t1 = [0, t75, t75, t75, t72; 0, -t73, -t73, -t73, -t74; 0, 0, 0, 0, -t84; 0, t74, t74, t74, t73; 0, t72, t72, t72, t75; 0, 0, 0, 0, -t83; 0, -t87, -t87, -t87, 0; 0, t76, t76, t76, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end