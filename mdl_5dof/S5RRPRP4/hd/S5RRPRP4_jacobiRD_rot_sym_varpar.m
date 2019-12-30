% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPRP4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_jacobiRD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:42:25
	% EndTime: 2019-12-29 18:42:25
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:42:25
	% EndTime: 2019-12-29 18:42:25
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
	% StartTime: 2019-12-29 18:42:25
	% EndTime: 2019-12-29 18:42:25
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
	% StartTime: 2019-12-29 18:42:20
	% EndTime: 2019-12-29 18:42:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (18->4), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t37 = qJ(1) + qJ(2);
	t36 = qJD(1) + qJD(2);
	t34 = t36 * cos(t37);
	t33 = t36 * sin(t37);
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; t34, t34, 0, 0, 0; t33, t33, 0, 0, 0; 0, 0, 0, 0, 0; -t33, -t33, 0, 0, 0; t34, t34, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:42:25
	% EndTime: 2019-12-29 18:42:25
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (61->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t70 = qJ(1) + qJ(2);
	t67 = sin(t70);
	t69 = qJD(1) + qJD(2);
	t78 = t69 * t67;
	t68 = cos(t70);
	t77 = t69 * t68;
	t71 = sin(qJ(4));
	t76 = t69 * t71;
	t72 = cos(qJ(4));
	t75 = t69 * t72;
	t74 = qJD(4) * t71;
	t73 = qJD(4) * t72;
	t66 = -t67 * t74 + t68 * t75;
	t65 = t67 * t73 + t68 * t76;
	t64 = t67 * t75 + t68 * t74;
	t63 = -t67 * t76 + t68 * t73;
	t1 = [t63, t63, 0, t66, 0; t65, t65, 0, t64, 0; 0, 0, 0, -t73, 0; -t64, -t64, 0, -t65, 0; t66, t66, 0, t63, 0; 0, 0, 0, t74, 0; -t77, -t77, 0, 0, 0; -t78, -t78, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:42:10
	% EndTime: 2019-12-29 18:42:10
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (62->16), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t223 = qJ(1) + qJ(2);
	t220 = sin(t223);
	t222 = qJD(1) + qJD(2);
	t231 = t222 * t220;
	t221 = cos(t223);
	t230 = t222 * t221;
	t224 = sin(qJ(4));
	t229 = t222 * t224;
	t225 = cos(qJ(4));
	t228 = t222 * t225;
	t227 = qJD(4) * t224;
	t226 = qJD(4) * t225;
	t217 = -t220 * t227 + t221 * t228;
	t216 = t220 * t226 + t221 * t229;
	t215 = t220 * t228 + t221 * t227;
	t214 = t220 * t229 - t221 * t226;
	t1 = [-t214, -t214, 0, t217, 0; t216, t216, 0, t215, 0; 0, 0, 0, -t226, 0; -t230, -t230, 0, 0, 0; -t231, -t231, 0, 0, 0; 0, 0, 0, 0, 0; t215, t215, 0, t216, 0; -t217, -t217, 0, t214, 0; 0, 0, 0, -t227, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end