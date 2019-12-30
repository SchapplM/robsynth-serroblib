% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRPR6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR6_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR6_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR6_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:02:20
	% EndTime: 2019-12-29 20:02:21
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:02:20
	% EndTime: 2019-12-29 20:02:21
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
	% StartTime: 2019-12-29 20:02:26
	% EndTime: 2019-12-29 20:02:26
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t34 = sin(qJ(1));
	t41 = qJD(1) * t34;
	t36 = cos(qJ(1));
	t40 = qJD(1) * t36;
	t33 = sin(qJ(2));
	t39 = qJD(2) * t33;
	t35 = cos(qJ(2));
	t38 = qJD(2) * t35;
	t37 = qJD(2) * t36;
	t32 = t34 * t39 - t35 * t40;
	t31 = t33 * t40 + t34 * t38;
	t30 = t33 * t37 + t35 * t41;
	t29 = t33 * t41 - t35 * t37;
	t1 = [t32, t29, 0, 0, 0; -t30, -t31, 0, 0, 0; 0, -t39, 0, 0, 0; t31, t30, 0, 0, 0; t29, t32, 0, 0, 0; 0, -t38, 0, 0, 0; -t41, 0, 0, 0, 0; t40, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:02:18
	% EndTime: 2019-12-29 20:02:18
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (61->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t71 = qJ(2) + qJ(3);
	t68 = sin(t71);
	t70 = qJD(2) + qJD(3);
	t79 = t70 * t68;
	t69 = cos(t71);
	t78 = t70 * t69;
	t72 = sin(qJ(1));
	t77 = t70 * t72;
	t73 = cos(qJ(1));
	t76 = t70 * t73;
	t75 = qJD(1) * t72;
	t74 = qJD(1) * t73;
	t67 = t68 * t77 - t69 * t74;
	t66 = t68 * t74 + t69 * t77;
	t65 = t68 * t76 + t69 * t75;
	t64 = t68 * t75 - t69 * t76;
	t1 = [t67, t64, t64, 0, 0; -t65, -t66, -t66, 0, 0; 0, -t79, -t79, 0, 0; t66, t65, t65, 0, 0; t64, t67, t67, 0, 0; 0, -t78, -t78, 0, 0; -t75, 0, 0, 0, 0; t74, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:02:23
	% EndTime: 2019-12-29 20:02:23
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (59->11), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t226 = qJ(2) + qJ(3);
	t223 = sin(t226);
	t225 = qJD(2) + qJD(3);
	t233 = t225 * t223;
	t227 = sin(qJ(1));
	t232 = t225 * t227;
	t228 = cos(qJ(1));
	t231 = t225 * t228;
	t230 = qJD(1) * t227;
	t229 = qJD(1) * t228;
	t224 = cos(t226);
	t222 = t225 * t224;
	t221 = -t223 * t232 + t224 * t229;
	t220 = -t223 * t229 - t224 * t232;
	t219 = -t223 * t231 - t224 * t230;
	t218 = t223 * t230 - t224 * t231;
	t1 = [-t221, t218, t218, 0, 0; t219, t220, t220, 0, 0; 0, -t233, -t233, 0, 0; -t230, 0, 0, 0, 0; t229, 0, 0, 0, 0; 0, 0, 0, 0, 0; t220, t219, t219, 0, 0; -t218, t221, t221, 0, 0; 0, t222, t222, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:02:21
	% EndTime: 2019-12-29 20:02:21
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (263->27), mult. (338->28), div. (0->0), fcn. (338->6), ass. (0->26)
	t167 = qJ(2) + qJ(3);
	t165 = cos(t167);
	t168 = sin(qJ(5));
	t184 = t165 * t168;
	t164 = sin(t167);
	t170 = cos(qJ(5));
	t185 = t164 * t170;
	t191 = t184 - t185;
	t166 = qJD(2) + qJD(3);
	t188 = qJD(5) - t166;
	t189 = -t164 * t168 - t165 * t170;
	t190 = t188 * t189;
	t187 = t191 * qJD(5) + t166 * t185;
	t171 = cos(qJ(1));
	t182 = t166 * t171;
	t169 = sin(qJ(1));
	t181 = qJD(1) * t169;
	t180 = qJD(1) * t171;
	t177 = t168 * t182;
	t172 = qJD(1) * t189;
	t155 = t188 * t191;
	t154 = -t171 * t172 + (t166 * t184 - t187) * t169;
	t153 = -t169 * t190 + t180 * t191;
	t152 = -t165 * t177 - t169 * t172 + t187 * t171;
	t151 = t164 * t177 - t181 * t185 + (t168 * t181 + t170 * t182) * t165 + t189 * qJD(5) * t171;
	t1 = [-t154, -t151, -t151, 0, t151; -t152, t153, t153, 0, -t153; 0, -t155, -t155, 0, t155; t153, -t152, -t152, 0, t152; t151, t154, t154, 0, -t154; 0, t190, t190, 0, -t190; t181, 0, 0, 0, 0; -t180, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end