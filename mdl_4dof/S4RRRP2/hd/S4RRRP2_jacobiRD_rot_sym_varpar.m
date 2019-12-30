% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% JRD_rot [9x4]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 13:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S4RRRP2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_jacobiRD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_jacobiRD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRRP2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_jacobiRD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:51:57
	% EndTime: 2019-12-29 13:51:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:51:49
	% EndTime: 2019-12-29 13:51:50
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0; -t31, 0, 0, 0; 0, 0, 0, 0; t31, 0, 0, 0; -t30, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:51:57
	% EndTime: 2019-12-29 13:51:57
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (22->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t47 = qJD(1) + qJD(2);
	t48 = qJ(1) + qJ(2);
	t49 = t47 * cos(t48);
	t44 = t47 * sin(t48);
	t1 = [-t49, -t49, 0, 0; -t44, -t44, 0, 0; 0, 0, 0, 0; t44, t44, 0, 0; -t49, -t49, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:51:52
	% EndTime: 2019-12-29 13:51:52
	% DurationCPUTime: 0.06s
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
	t1 = [t65, t65, t62, 0; -t63, -t63, -t64, 0; 0, 0, -t74, 0; t64, t64, t63, 0; t62, t62, t65, 0; 0, 0, -t73, 0; -t77, -t77, 0, 0; t66, t66, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:51:45
	% EndTime: 2019-12-29 13:51:45
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (60->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t81 = qJ(1) + qJ(2);
	t78 = sin(t81);
	t80 = qJD(1) + qJD(2);
	t88 = t80 * t78;
	t82 = sin(qJ(3));
	t87 = t80 * t82;
	t83 = cos(qJ(3));
	t86 = t80 * t83;
	t85 = qJD(3) * t82;
	t84 = qJD(3) * t83;
	t79 = cos(t81);
	t77 = t80 * t79;
	t76 = t78 * t85 - t79 * t86;
	t75 = t78 * t84 + t79 * t87;
	t74 = t78 * t86 + t79 * t85;
	t73 = t78 * t87 - t79 * t84;
	t1 = [t76, t76, t73, 0; -t74, -t74, -t75, 0; 0, 0, -t85, 0; t75, t75, t74, 0; t73, t73, t76, 0; 0, 0, -t84, 0; -t88, -t88, 0, 0; t77, t77, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,4);
end