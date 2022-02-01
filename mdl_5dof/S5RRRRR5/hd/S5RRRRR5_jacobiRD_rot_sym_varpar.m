% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR5
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
%   Siehe auch: S5RRRRR5_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRRR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:02:45
	% EndTime: 2022-01-20 12:02:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:02:45
	% EndTime: 2022-01-20 12:02:46
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:02:46
	% EndTime: 2022-01-20 12:02:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t47 = qJD(1) + qJD(2);
	t48 = qJ(1) + qJ(2);
	t49 = t47 * cos(t48);
	t44 = t47 * sin(t48);
	t1 = [-t49, -t49, 0, 0, 0; -t44, -t44, 0, 0, 0; 0, 0, 0, 0, 0; t44, t44, 0, 0, 0; -t49, -t49, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:02:45
	% EndTime: 2022-01-20 12:02:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (57->11), mult. (12->2), div. (0->0), fcn. (12->2), ass. (0->5)
	t59 = qJD(1) + qJD(2) + qJD(3);
	t60 = qJ(1) + qJ(2) + qJ(3);
	t61 = t59 * cos(t60);
	t56 = t59 * sin(t60);
	t1 = [-t61, -t61, -t61, 0, 0; -t56, -t56, -t56, 0, 0; 0, 0, 0, 0, 0; t56, t56, t56, 0, 0; -t61, -t61, -t61, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:02:46
	% EndTime: 2022-01-20 12:02:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (141->15), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t82 = qJ(1) + qJ(2) + qJ(3);
	t79 = sin(t82);
	t81 = qJD(1) + qJD(2) + qJD(3);
	t89 = t81 * t79;
	t83 = sin(qJ(4));
	t88 = t81 * t83;
	t84 = cos(qJ(4));
	t87 = t81 * t84;
	t86 = qJD(4) * t83;
	t85 = qJD(4) * t84;
	t80 = cos(t82);
	t78 = t81 * t80;
	t77 = t79 * t86 - t80 * t87;
	t76 = t79 * t85 + t80 * t88;
	t75 = t79 * t87 + t80 * t86;
	t74 = t79 * t88 - t80 * t85;
	t1 = [t77, t77, t77, t74, 0; -t75, -t75, -t75, -t76, 0; 0, 0, 0, -t86, 0; t76, t76, t76, t75, 0; t74, t74, t74, t77, 0; 0, 0, 0, -t85, 0; -t89, -t89, -t89, 0, 0; t78, t78, t78, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:02:46
	% EndTime: 2022-01-20 12:02:46
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (239->20), mult. (90->14), div. (0->0), fcn. (90->4), ass. (0->19)
	t113 = qJ(1) + qJ(2) + qJ(3);
	t108 = sin(t113);
	t110 = qJD(1) + qJD(2) + qJD(3);
	t120 = t110 * t108;
	t115 = qJ(4) + qJ(5);
	t111 = sin(t115);
	t119 = t110 * t111;
	t112 = cos(t115);
	t118 = t110 * t112;
	t114 = qJD(4) + qJD(5);
	t117 = t114 * t111;
	t116 = t114 * t112;
	t109 = cos(t113);
	t107 = t110 * t109;
	t106 = t108 * t117 - t109 * t118;
	t105 = t108 * t116 + t109 * t119;
	t104 = t108 * t118 + t109 * t117;
	t103 = t108 * t119 - t109 * t116;
	t1 = [t106, t106, t106, t103, t103; -t104, -t104, -t104, -t105, -t105; 0, 0, 0, -t117, -t117; t105, t105, t105, t104, t104; t103, t103, t103, t106, t106; 0, 0, 0, -t116, -t116; -t120, -t120, -t120, 0, 0; t107, t107, t107, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end