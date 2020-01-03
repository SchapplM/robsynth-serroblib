% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPRP2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:50:06
	% EndTime: 2019-12-31 19:50:06
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:50:06
	% EndTime: 2019-12-31 19:50:06
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
	% StartTime: 2019-12-31 19:50:06
	% EndTime: 2019-12-31 19:50:06
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
	% StartTime: 2019-12-31 19:50:06
	% EndTime: 2019-12-31 19:50:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t53 = qJ(1) + qJ(2) + pkin(8);
	t54 = qJD(1) + qJD(2);
	t55 = t54 * cos(t53);
	t50 = t54 * sin(t53);
	t1 = [-t55, -t55, 0, 0, 0; -t50, -t50, 0, 0, 0; 0, 0, 0, 0, 0; t50, t50, 0, 0, 0; -t55, -t55, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:50:06
	% EndTime: 2019-12-31 19:50:06
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (88->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t74 = qJ(1) + qJ(2) + pkin(8);
	t72 = sin(t74);
	t75 = qJD(1) + qJD(2);
	t82 = t75 * t72;
	t76 = sin(qJ(4));
	t81 = t75 * t76;
	t77 = cos(qJ(4));
	t80 = t75 * t77;
	t79 = qJD(4) * t76;
	t78 = qJD(4) * t77;
	t73 = cos(t74);
	t71 = t75 * t73;
	t70 = t72 * t79 - t73 * t80;
	t69 = t72 * t78 + t73 * t81;
	t68 = t72 * t80 + t73 * t79;
	t67 = t72 * t81 - t73 * t78;
	t1 = [t70, t70, 0, t67, 0; -t68, -t68, 0, -t69, 0; 0, 0, 0, -t79, 0; t69, t69, 0, t68, 0; t67, t67, 0, t70, 0; 0, 0, 0, -t78, 0; -t82, -t82, 0, 0, 0; t71, t71, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:50:07
	% EndTime: 2019-12-31 19:50:07
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (87->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t209 = qJ(1) + qJ(2) + pkin(8);
	t207 = sin(t209);
	t210 = qJD(1) + qJD(2);
	t217 = t210 * t207;
	t211 = sin(qJ(4));
	t216 = t210 * t211;
	t212 = cos(qJ(4));
	t215 = t210 * t212;
	t214 = qJD(4) * t211;
	t213 = qJD(4) * t212;
	t208 = cos(t209);
	t206 = t210 * t208;
	t203 = -t207 * t214 + t208 * t215;
	t202 = -t207 * t213 - t208 * t216;
	t201 = -t207 * t215 - t208 * t214;
	t200 = t207 * t216 - t208 * t213;
	t1 = [-t203, -t203, 0, t200, 0; t201, t201, 0, t202, 0; 0, 0, 0, -t214, 0; -t217, -t217, 0, 0, 0; t206, t206, 0, 0, 0; 0, 0, 0, 0, 0; t202, t202, 0, t201, 0; -t200, -t200, 0, t203, 0; 0, 0, 0, t213, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end