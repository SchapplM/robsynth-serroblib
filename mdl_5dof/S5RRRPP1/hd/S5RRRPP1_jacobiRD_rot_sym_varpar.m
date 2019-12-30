% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRPP1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:33:47
	% EndTime: 2019-12-29 19:33:47
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:33:37
	% EndTime: 2019-12-29 19:33:37
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:33:42
	% EndTime: 2019-12-29 19:33:42
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-12-29 19:33:32
	% EndTime: 2019-12-29 19:33:32
	% DurationCPUTime: 0.07s
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
	% StartTime: 2019-12-29 19:33:42
	% EndTime: 2019-12-29 19:33:42
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (86->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->16)
	t90 = qJ(1) + qJ(2);
	t86 = sin(t90);
	t88 = qJD(1) + qJD(2);
	t93 = t88 * t86;
	t87 = cos(t90);
	t83 = t88 * t87;
	t92 = qJD(3) * t86;
	t91 = qJD(3) * t87;
	t89 = qJ(3) + pkin(8);
	t85 = cos(t89);
	t84 = sin(t89);
	t82 = -t85 * t83 + t84 * t92;
	t81 = t84 * t83 + t85 * t92;
	t80 = t84 * t91 + t85 * t93;
	t79 = t84 * t93 - t85 * t91;
	t1 = [t82, t82, t79, 0, 0; -t80, -t80, -t81, 0, 0; 0, 0, -qJD(3) * t84, 0, 0; t81, t81, t80, 0, 0; t79, t79, t82, 0, 0; 0, 0, -qJD(3) * t85, 0, 0; -t93, -t93, 0, 0, 0; t83, t83, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:33:44
	% EndTime: 2019-12-29 19:33:44
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (85->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->16)
	t222 = qJ(1) + qJ(2);
	t218 = sin(t222);
	t220 = qJD(1) + qJD(2);
	t225 = t220 * t218;
	t219 = cos(t222);
	t215 = t220 * t219;
	t224 = qJD(3) * t218;
	t223 = qJD(3) * t219;
	t221 = qJ(3) + pkin(8);
	t217 = cos(t221);
	t216 = sin(t221);
	t212 = t217 * t215 - t216 * t224;
	t211 = -t216 * t215 - t217 * t224;
	t210 = -t216 * t223 - t217 * t225;
	t209 = t216 * t225 - t217 * t223;
	t1 = [-t212, -t212, t209, 0, 0; t210, t210, t211, 0, 0; 0, 0, -qJD(3) * t216, 0, 0; -t225, -t225, 0, 0, 0; t215, t215, 0, 0, 0; 0, 0, 0, 0, 0; t211, t211, t210, 0, 0; -t209, -t209, t212, 0, 0; 0, 0, qJD(3) * t217, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end