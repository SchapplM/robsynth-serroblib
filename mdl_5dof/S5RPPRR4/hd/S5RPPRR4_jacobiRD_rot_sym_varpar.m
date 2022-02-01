% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRR4
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
%   Siehe auch: S5RPPRR4_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPRR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
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
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t22 = qJD(1) * sin(qJ(1));
	t21 = qJD(1) * cos(qJ(1));
	t18 = cos(pkin(8));
	t17 = sin(pkin(8));
	t1 = [-t18 * t21, 0, 0, 0, 0; -t18 * t22, 0, 0, 0, 0; 0, 0, 0, 0, 0; t17 * t21, 0, 0, 0, 0; t17 * t22, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t22, 0, 0, 0, 0; t21, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:42
	% EndTime: 2022-01-23 09:17:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t128 = cos(pkin(8));
	t129 = sin(qJ(1));
	t133 = t128 * t129;
	t130 = cos(qJ(1));
	t132 = t128 * t130;
	t131 = qJD(1) * sin(pkin(8));
	t127 = cos(pkin(9));
	t125 = sin(pkin(9));
	t1 = [(-t125 * t129 - t127 * t132) * qJD(1), 0, 0, 0, 0; (t125 * t130 - t127 * t133) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; (t125 * t132 - t127 * t129) * qJD(1), 0, 0, 0, 0; (t125 * t133 + t127 * t130) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; -t130 * t131, 0, 0, 0, 0; -t129 * t131, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:42
	% EndTime: 2022-01-23 09:17:42
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (61->14), mult. (88->24), div. (0->0), fcn. (88->6), ass. (0->20)
	t169 = cos(pkin(8));
	t170 = sin(qJ(1));
	t179 = t169 * t170;
	t171 = cos(qJ(1));
	t178 = t169 * t171;
	t168 = sin(pkin(8));
	t177 = qJD(1) * t168;
	t176 = qJD(4) * t168;
	t167 = pkin(9) + qJ(4);
	t165 = sin(t167);
	t166 = cos(t167);
	t175 = -t165 * t170 - t166 * t178;
	t174 = -t165 * t171 + t166 * t179;
	t173 = t165 * t178 - t166 * t170;
	t172 = t165 * t179 + t166 * t171;
	t164 = t175 * qJD(1) + t172 * qJD(4);
	t163 = t173 * qJD(1) + t174 * qJD(4);
	t162 = t174 * qJD(1) + t173 * qJD(4);
	t161 = t172 * qJD(1) + t175 * qJD(4);
	t1 = [t164, 0, 0, t161, 0; -t162, 0, 0, -t163, 0; 0, 0, 0, -t166 * t176, 0; t163, 0, 0, t162, 0; t161, 0, 0, t164, 0; 0, 0, 0, t165 * t176, 0; -t171 * t177, 0, 0, 0, 0; -t170 * t177, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:42
	% EndTime: 2022-01-23 09:17:42
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (172->19), mult. (132->28), div. (0->0), fcn. (132->6), ass. (0->24)
	t213 = pkin(9) + qJ(4) + qJ(5);
	t211 = sin(t213);
	t218 = cos(qJ(1));
	t228 = t211 * t218;
	t212 = cos(t213);
	t217 = sin(qJ(1));
	t227 = t212 * t217;
	t214 = qJD(4) + qJD(5);
	t215 = sin(pkin(8));
	t226 = t214 * t215;
	t216 = cos(pkin(8));
	t225 = t216 * t217;
	t224 = t216 * t218;
	t223 = qJD(1) * t217;
	t222 = qJD(1) * t218;
	t221 = t212 * t226;
	t220 = -t211 * t217 - t212 * t224;
	t219 = t211 * t225 + t212 * t218;
	t208 = t211 * t226;
	t207 = t220 * qJD(1) + t219 * t214;
	t206 = -t214 * t228 - t212 * t223 + (t211 * t222 + t214 * t227) * t216;
	t205 = (t211 * t224 - t227) * t214 + (t212 * t225 - t228) * qJD(1);
	t204 = t219 * qJD(1) + t220 * t214;
	t1 = [t207, 0, 0, t204, t204; -t205, 0, 0, -t206, -t206; 0, 0, 0, -t221, -t221; t206, 0, 0, t205, t205; t204, 0, 0, t207, t207; 0, 0, 0, t208, t208; -t215 * t222, 0, 0, 0, 0; -t215 * t223, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end