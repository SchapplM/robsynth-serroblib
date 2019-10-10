% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:52
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPR13_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR13_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:52:10
	% EndTime: 2019-10-10 12:52:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:52:10
	% EndTime: 2019-10-10 12:52:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:52:10
	% EndTime: 2019-10-10 12:52:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:52:10
	% EndTime: 2019-10-10 12:52:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t103 = sin(qJ(2));
	t104 = sin(qJ(1));
	t111 = t103 * t104;
	t106 = cos(qJ(1));
	t110 = t103 * t106;
	t105 = cos(qJ(2));
	t109 = t104 * t105;
	t108 = t105 * t106;
	t101 = sin(pkin(6));
	t107 = qJD(1) * t101;
	t102 = cos(pkin(6));
	t1 = [0, t106 * t107, (-t102 * t111 + t108) * qJD(2) + (t102 * t108 - t111) * qJD(1), 0, 0, 0; 0, t104 * t107, (t102 * t110 + t109) * qJD(2) + (t102 * t109 + t110) * qJD(1), 0, 0, 0; 0, 0, t101 * qJD(2) * t103, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:52:10
	% EndTime: 2019-10-10 12:52:10
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (22->16), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->21)
	t155 = sin(pkin(6));
	t160 = cos(qJ(3));
	t174 = t155 * t160;
	t162 = cos(qJ(1));
	t173 = t155 * t162;
	t158 = sin(qJ(2));
	t159 = sin(qJ(1));
	t172 = t158 * t159;
	t171 = t158 * t162;
	t161 = cos(qJ(2));
	t170 = t159 * t161;
	t169 = t162 * t161;
	t168 = qJD(1) * t155;
	t157 = sin(qJ(3));
	t167 = qJD(2) * t157;
	t156 = cos(pkin(6));
	t166 = t156 * t169 - t172;
	t165 = t156 * t170 + t171;
	t164 = t156 * t171 + t170;
	t163 = -t156 * t172 + t169;
	t1 = [0, t162 * t168, t166 * qJD(1) + t163 * qJD(2), (t159 * t155 * t157 + t163 * t160) * qJD(3) - t165 * t167 + (-t164 * t157 - t160 * t173) * qJD(1), 0, 0; 0, t159 * t168, t165 * qJD(1) + t164 * qJD(2), (-t157 * t173 + t164 * t160) * qJD(3) + t166 * t167 + (t163 * t157 - t159 * t174) * qJD(1), 0, 0; 0, 0, t155 * qJD(2) * t158, t155 * t161 * t167 + (t156 * t157 + t158 * t174) * qJD(3), 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:52:11
	% EndTime: 2019-10-10 12:52:11
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (22->16), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->21)
	t188 = sin(pkin(6));
	t193 = cos(qJ(3));
	t207 = t188 * t193;
	t195 = cos(qJ(1));
	t206 = t188 * t195;
	t191 = sin(qJ(2));
	t192 = sin(qJ(1));
	t205 = t191 * t192;
	t204 = t191 * t195;
	t194 = cos(qJ(2));
	t203 = t192 * t194;
	t202 = t195 * t194;
	t201 = qJD(1) * t188;
	t190 = sin(qJ(3));
	t200 = qJD(2) * t190;
	t189 = cos(pkin(6));
	t199 = t189 * t202 - t205;
	t198 = t189 * t203 + t204;
	t197 = t189 * t204 + t203;
	t196 = -t189 * t205 + t202;
	t1 = [0, t195 * t201, t199 * qJD(1) + t196 * qJD(2), (t192 * t188 * t190 + t196 * t193) * qJD(3) - t198 * t200 + (-t197 * t190 - t193 * t206) * qJD(1), 0, 0; 0, t192 * t201, t198 * qJD(1) + t197 * qJD(2), (-t190 * t206 + t197 * t193) * qJD(3) + t199 * t200 + (t196 * t190 - t192 * t207) * qJD(1), 0, 0; 0, 0, t188 * qJD(2) * t191, t188 * t194 * t200 + (t189 * t190 + t191 * t207) * qJD(3), 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:52:11
	% EndTime: 2019-10-10 12:52:11
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (38->19), mult. (130->40), div. (0->0), fcn. (134->8), ass. (0->24)
	t217 = sin(pkin(6));
	t222 = cos(qJ(3));
	t236 = t217 * t222;
	t224 = cos(qJ(1));
	t235 = t217 * t224;
	t220 = sin(qJ(2));
	t221 = sin(qJ(1));
	t234 = t220 * t221;
	t233 = t220 * t224;
	t223 = cos(qJ(2));
	t232 = t221 * t223;
	t231 = t224 * t223;
	t230 = qJD(1) * t217;
	t219 = sin(qJ(3));
	t229 = qJD(2) * t219;
	t218 = cos(pkin(6));
	t228 = t218 * t231 - t234;
	t227 = t218 * t232 + t233;
	t226 = t218 * t233 + t232;
	t225 = -t218 * t234 + t231;
	t216 = t217 * t223 * t229 + (t218 * t219 + t220 * t236) * qJD(3);
	t215 = (-t219 * t235 + t226 * t222) * qJD(3) + t228 * t229 + (t225 * t219 - t221 * t236) * qJD(1);
	t214 = (t221 * t217 * t219 + t225 * t222) * qJD(3) - t227 * t229 + (-t226 * t219 - t222 * t235) * qJD(1);
	t1 = [0, t224 * t230, t228 * qJD(1) + t225 * qJD(2), t214, 0, -t214; 0, t221 * t230, t227 * qJD(1) + t226 * qJD(2), t215, 0, -t215; 0, 0, t217 * qJD(2) * t220, t216, 0, -t216;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end