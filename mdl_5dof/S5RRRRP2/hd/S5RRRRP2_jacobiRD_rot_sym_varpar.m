% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRP2
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
%   Siehe auch: S5RRRRP2_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRRP2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:49:58
	% EndTime: 2022-01-20 11:49:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:49:58
	% EndTime: 2022-01-20 11:49:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:49:58
	% EndTime: 2022-01-20 11:49:58
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
	% StartTime: 2022-01-20 11:49:58
	% EndTime: 2022-01-20 11:49:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (60->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t70 = qJ(1) + qJ(2);
	t67 = sin(t70);
	t69 = qJD(1) + qJD(2);
	t77 = t69 * t67;
	t71 = sin(qJ(3));
	t76 = t69 * t71;
	t72 = cos(qJ(3));
	t75 = t69 * t72;
	t74 = qJD(3) * t71;
	t73 = qJD(3) * t72;
	t68 = cos(t70);
	t66 = t69 * t68;
	t65 = t67 * t74 - t68 * t75;
	t64 = t67 * t73 + t68 * t76;
	t63 = t67 * t75 + t68 * t74;
	t62 = t67 * t76 - t68 * t73;
	t1 = [t65, t65, t62, 0, 0; -t63, -t63, -t64, 0, 0; 0, 0, -t74, 0, 0; t64, t64, t63, 0, 0; t62, t62, t65, 0, 0; 0, 0, -t73, 0, 0; -t77, -t77, 0, 0, 0; t66, t66, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:49:58
	% EndTime: 2022-01-20 11:49:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (134->18), mult. (72->12), div. (0->0), fcn. (72->4), ass. (0->17)
	t101 = qJ(3) + qJ(4);
	t95 = sin(t101);
	t99 = qJD(3) + qJD(4);
	t105 = t99 * t95;
	t97 = cos(t101);
	t104 = t99 * t97;
	t100 = qJD(1) + qJD(2);
	t102 = qJ(1) + qJ(2);
	t96 = sin(t102);
	t103 = t100 * t96;
	t98 = cos(t102);
	t94 = t100 * t98;
	t93 = t96 * t105 - t97 * t94;
	t92 = t96 * t104 + t95 * t94;
	t91 = t97 * t103 + t98 * t105;
	t90 = t95 * t103 - t98 * t104;
	t1 = [t93, t93, t90, t90, 0; -t91, -t91, -t92, -t92, 0; 0, 0, -t105, -t105, 0; t92, t92, t91, t91, 0; t90, t90, t93, t93, 0; 0, 0, -t104, -t104, 0; -t103, -t103, 0, 0, 0; t94, t94, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:49:58
	% EndTime: 2022-01-20 11:49:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (134->18), mult. (72->12), div. (0->0), fcn. (72->4), ass. (0->17)
	t111 = qJ(3) + qJ(4);
	t105 = sin(t111);
	t109 = qJD(3) + qJD(4);
	t115 = t109 * t105;
	t107 = cos(t111);
	t114 = t109 * t107;
	t112 = qJ(1) + qJ(2);
	t106 = sin(t112);
	t110 = qJD(1) + qJD(2);
	t113 = t110 * t106;
	t108 = cos(t112);
	t104 = t110 * t108;
	t103 = -t107 * t104 + t106 * t115;
	t102 = t105 * t104 + t106 * t114;
	t101 = t107 * t113 + t108 * t115;
	t100 = t105 * t113 - t108 * t114;
	t1 = [t103, t103, t100, t100, 0; -t101, -t101, -t102, -t102, 0; 0, 0, -t115, -t115, 0; t102, t102, t101, t101, 0; t100, t100, t103, t103, 0; 0, 0, -t114, -t114, 0; -t113, -t113, 0, 0, 0; t104, t104, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end