% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRP4
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
%   Siehe auch: S5RPRRP4_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRRP4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
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
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
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
	% StartTime: 2022-01-23 09:33:25
	% EndTime: 2022-01-23 09:33:25
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
	% StartTime: 2022-01-23 09:33:25
	% EndTime: 2022-01-23 09:33:25
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (120->19), mult. (132->28), div. (0->0), fcn. (132->6), ass. (0->24)
	t209 = qJ(3) + qJ(4);
	t206 = sin(t209);
	t213 = cos(qJ(1));
	t223 = t206 * t213;
	t207 = cos(t209);
	t212 = sin(qJ(1));
	t222 = t207 * t212;
	t208 = qJD(3) + qJD(4);
	t210 = sin(pkin(8));
	t221 = t208 * t210;
	t211 = cos(pkin(8));
	t220 = t211 * t212;
	t219 = t211 * t213;
	t218 = qJD(1) * t212;
	t217 = qJD(1) * t213;
	t216 = t207 * t221;
	t215 = -t206 * t212 - t207 * t219;
	t214 = t206 * t220 + t207 * t213;
	t203 = t206 * t221;
	t202 = t215 * qJD(1) + t214 * t208;
	t201 = -t208 * t223 - t207 * t218 + (t206 * t217 + t208 * t222) * t211;
	t200 = (t206 * t219 - t222) * t208 + (t207 * t220 - t223) * qJD(1);
	t199 = t214 * qJD(1) + t215 * t208;
	t1 = [t202, 0, t199, t199, 0; -t200, 0, -t201, -t201, 0; 0, 0, -t216, -t216, 0; t201, 0, t200, t200, 0; t199, 0, t202, t202, 0; 0, 0, t203, t203, 0; -t210 * t217, 0, 0, 0, 0; -t210 * t218, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:25
	% EndTime: 2022-01-23 09:33:25
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (120->19), mult. (132->28), div. (0->0), fcn. (132->6), ass. (0->24)
	t218 = qJ(3) + qJ(4);
	t215 = sin(t218);
	t222 = cos(qJ(1));
	t232 = t215 * t222;
	t216 = cos(t218);
	t221 = sin(qJ(1));
	t231 = t216 * t221;
	t217 = qJD(3) + qJD(4);
	t219 = sin(pkin(8));
	t230 = t217 * t219;
	t220 = cos(pkin(8));
	t229 = t220 * t221;
	t228 = t220 * t222;
	t227 = qJD(1) * t221;
	t226 = qJD(1) * t222;
	t225 = t216 * t230;
	t224 = -t215 * t221 - t216 * t228;
	t223 = t215 * t229 + t216 * t222;
	t212 = t215 * t230;
	t211 = t224 * qJD(1) + t223 * t217;
	t210 = -t217 * t232 - t216 * t227 + (t215 * t226 + t217 * t231) * t220;
	t209 = (t215 * t228 - t231) * t217 + (t216 * t229 - t232) * qJD(1);
	t208 = t223 * qJD(1) + t224 * t217;
	t1 = [t211, 0, t208, t208, 0; -t209, 0, -t210, -t210, 0; 0, 0, -t225, -t225, 0; t210, 0, t209, t209, 0; t208, 0, t211, t211, 0; 0, 0, t212, t212, 0; -t219 * t226, 0, 0, 0, 0; -t219 * t227, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end