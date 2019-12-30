% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRP5
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
% Datum: 2019-12-29 17:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRRP5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:16:21
	% EndTime: 2019-12-29 17:16:21
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:16:21
	% EndTime: 2019-12-29 17:16:22
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
	% StartTime: 2019-12-29 17:16:17
	% EndTime: 2019-12-29 17:16:17
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(8);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0; -t37, 0, 0, 0, 0; 0, 0, 0, 0, 0; t37, 0, 0, 0, 0; -t36, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:16:16
	% EndTime: 2019-12-29 17:16:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (30->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t48 = qJ(1) + pkin(8) + qJ(3);
	t49 = qJD(1) + qJD(3);
	t50 = t49 * cos(t48);
	t45 = t49 * sin(t48);
	t1 = [-t50, 0, -t50, 0, 0; -t45, 0, -t45, 0, 0; 0, 0, 0, 0, 0; t45, 0, t45, 0, 0; -t50, 0, -t50, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:16:22
	% EndTime: 2019-12-29 17:16:22
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (88->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t70 = qJ(1) + pkin(8) + qJ(3);
	t68 = sin(t70);
	t71 = qJD(1) + qJD(3);
	t78 = t71 * t68;
	t72 = sin(qJ(4));
	t77 = t71 * t72;
	t73 = cos(qJ(4));
	t76 = t71 * t73;
	t75 = qJD(4) * t72;
	t74 = qJD(4) * t73;
	t69 = cos(t70);
	t67 = t71 * t69;
	t66 = t68 * t75 - t69 * t76;
	t65 = t68 * t74 + t69 * t77;
	t64 = t68 * t76 + t69 * t75;
	t63 = t68 * t77 - t69 * t74;
	t1 = [t66, 0, t66, t63, 0; -t64, 0, -t64, -t65, 0; 0, 0, 0, -t75, 0; t65, 0, t65, t64, 0; t63, 0, t63, t66, 0; 0, 0, 0, -t74, 0; -t78, 0, -t78, 0, 0; t67, 0, t67, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:16:23
	% EndTime: 2019-12-29 17:16:23
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (87->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t206 = qJ(1) + pkin(8) + qJ(3);
	t204 = sin(t206);
	t207 = qJD(1) + qJD(3);
	t214 = t207 * t204;
	t208 = sin(qJ(4));
	t213 = t207 * t208;
	t209 = cos(qJ(4));
	t212 = t207 * t209;
	t211 = qJD(4) * t208;
	t210 = qJD(4) * t209;
	t205 = cos(t206);
	t203 = t207 * t205;
	t200 = -t204 * t211 + t205 * t212;
	t199 = -t204 * t210 - t205 * t213;
	t198 = -t204 * t212 - t205 * t211;
	t197 = t204 * t213 - t205 * t210;
	t1 = [-t200, 0, -t200, t197, 0; t198, 0, t198, t199, 0; 0, 0, 0, -t211, 0; -t214, 0, -t214, 0, 0; t203, 0, t203, 0, 0; 0, 0, 0, 0, 0; t199, 0, t199, t198, 0; -t197, 0, -t197, t200, 0; 0, 0, 0, t210, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end