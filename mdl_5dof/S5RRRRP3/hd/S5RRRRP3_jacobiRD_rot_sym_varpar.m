% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRRP3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
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
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
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
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.07s
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
	% StartTime: 2019-12-29 20:25:02
	% EndTime: 2019-12-29 20:25:02
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (140->16), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t217 = qJ(1) + qJ(2) + qJ(3);
	t214 = sin(t217);
	t216 = qJD(1) + qJD(2) + qJD(3);
	t224 = t216 * t214;
	t218 = sin(qJ(4));
	t223 = t216 * t218;
	t219 = cos(qJ(4));
	t222 = t216 * t219;
	t221 = qJD(4) * t218;
	t220 = qJD(4) * t219;
	t215 = cos(t217);
	t211 = t216 * t215;
	t210 = -t214 * t221 + t215 * t222;
	t209 = -t214 * t220 - t215 * t223;
	t208 = -t214 * t222 - t215 * t221;
	t207 = t214 * t223 - t215 * t220;
	t1 = [-t210, -t210, -t210, t207, 0; t208, t208, t208, t209, 0; 0, 0, 0, -t221, 0; -t224, -t224, -t224, 0, 0; t211, t211, t211, 0, 0; 0, 0, 0, 0, 0; t209, t209, t209, t208, 0; -t207, -t207, -t207, t210, 0; 0, 0, 0, t220, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end