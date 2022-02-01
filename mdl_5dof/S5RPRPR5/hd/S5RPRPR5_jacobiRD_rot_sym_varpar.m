% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPR5
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
%   Siehe auch: S5RPRPR5_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRPR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
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
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
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
	% StartTime: 2022-01-23 09:26:47
	% EndTime: 2022-01-23 09:26:47
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (27->13), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t162 = sin(qJ(3));
	t163 = sin(qJ(1));
	t175 = t162 * t163;
	t165 = cos(qJ(1));
	t174 = t162 * t165;
	t164 = cos(qJ(3));
	t173 = t163 * t164;
	t172 = t164 * t165;
	t160 = sin(pkin(8));
	t171 = qJD(1) * t160;
	t170 = qJD(3) * t160;
	t161 = cos(pkin(8));
	t169 = -t161 * t172 - t175;
	t168 = t161 * t173 - t174;
	t167 = t161 * t174 - t173;
	t166 = t161 * t175 + t172;
	t159 = t169 * qJD(1) + t166 * qJD(3);
	t158 = t167 * qJD(1) + t168 * qJD(3);
	t157 = t168 * qJD(1) + t167 * qJD(3);
	t156 = t166 * qJD(1) + t169 * qJD(3);
	t1 = [t159, 0, t156, 0, 0; -t157, 0, -t158, 0, 0; 0, 0, -t164 * t170, 0, 0; t158, 0, t157, 0, 0; t156, 0, t159, 0, 0; 0, 0, t162 * t170, 0, 0; -t165 * t171, 0, 0, 0, 0; -t163 * t171, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:47
	% EndTime: 2022-01-23 09:26:47
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (61->14), mult. (88->24), div. (0->0), fcn. (88->6), ass. (0->20)
	t182 = cos(pkin(8));
	t183 = sin(qJ(1));
	t192 = t182 * t183;
	t184 = cos(qJ(1));
	t191 = t182 * t184;
	t181 = sin(pkin(8));
	t190 = qJD(1) * t181;
	t189 = qJD(3) * t181;
	t180 = qJ(3) + pkin(9);
	t178 = sin(t180);
	t179 = cos(t180);
	t188 = -t178 * t183 - t179 * t191;
	t187 = -t178 * t184 + t179 * t192;
	t186 = t178 * t191 - t179 * t183;
	t185 = t178 * t192 + t179 * t184;
	t177 = t188 * qJD(1) + t185 * qJD(3);
	t176 = t186 * qJD(1) + t187 * qJD(3);
	t175 = t187 * qJD(1) + t186 * qJD(3);
	t174 = t185 * qJD(1) + t188 * qJD(3);
	t1 = [t177, 0, t174, 0, 0; -t175, 0, -t176, 0, 0; 0, 0, -t179 * t189, 0, 0; t176, 0, t175, 0, 0; t174, 0, t177, 0, 0; 0, 0, t178 * t189, 0, 0; -t184 * t190, 0, 0, 0, 0; -t183 * t190, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:47
	% EndTime: 2022-01-23 09:26:47
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (172->19), mult. (132->28), div. (0->0), fcn. (132->6), ass. (0->24)
	t218 = qJ(3) + pkin(9) + qJ(5);
	t216 = sin(t218);
	t223 = cos(qJ(1));
	t233 = t216 * t223;
	t217 = cos(t218);
	t222 = sin(qJ(1));
	t232 = t217 * t222;
	t219 = qJD(3) + qJD(5);
	t220 = sin(pkin(8));
	t231 = t219 * t220;
	t221 = cos(pkin(8));
	t230 = t221 * t222;
	t229 = t221 * t223;
	t228 = qJD(1) * t222;
	t227 = qJD(1) * t223;
	t226 = t217 * t231;
	t225 = -t216 * t222 - t217 * t229;
	t224 = t216 * t230 + t217 * t223;
	t213 = t216 * t231;
	t212 = t225 * qJD(1) + t224 * t219;
	t211 = -t219 * t233 - t217 * t228 + (t216 * t227 + t219 * t232) * t221;
	t210 = (t216 * t229 - t232) * t219 + (t217 * t230 - t233) * qJD(1);
	t209 = t224 * qJD(1) + t225 * t219;
	t1 = [t212, 0, t209, 0, t209; -t210, 0, -t211, 0, -t211; 0, 0, -t226, 0, -t226; t211, 0, t210, 0, t210; t209, 0, t212, 0, t212; 0, 0, t213, 0, t213; -t220 * t227, 0, 0, 0, 0; -t220 * t228, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end