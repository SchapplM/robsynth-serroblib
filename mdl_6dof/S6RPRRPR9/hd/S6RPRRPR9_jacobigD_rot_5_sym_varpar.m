% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRPR9_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:29
% EndTime: 2019-02-26 21:05:29
% DurationCPUTime: 0.06s
% Computational Cost: add. (24->18), mult. (92->42), div. (0->0), fcn. (96->10), ass. (0->26)
t196 = sin(pkin(7));
t197 = sin(pkin(6));
t217 = t197 * t196;
t199 = cos(pkin(7));
t203 = cos(qJ(3));
t216 = t199 * t203;
t195 = sin(pkin(12));
t202 = sin(qJ(1));
t215 = t202 * t195;
t198 = cos(pkin(12));
t214 = t202 * t198;
t204 = cos(qJ(1));
t213 = t204 * t195;
t212 = t204 * t198;
t211 = t202 * t217;
t210 = t204 * t217;
t209 = qJD(1) * t197 * t199;
t200 = cos(pkin(6));
t208 = t200 * t212 - t215;
t207 = -t200 * t214 - t213;
t206 = t200 * t213 + t214;
t205 = -t200 * t215 + t212;
t201 = sin(qJ(3));
t194 = t207 * qJD(1);
t193 = t208 * qJD(1);
t1 = [0, 0, t193 * t196 + t204 * t209, t193 * t216 + (t205 * t203 + (t207 * t199 + t211) * t201) * qJD(3) + (-t206 * t201 - t203 * t210) * qJD(1), 0, 0; 0, 0, -t194 * t196 + t202 * t209, -t194 * t216 + (t206 * t203 + (t208 * t199 - t210) * t201) * qJD(3) + (t205 * t201 - t203 * t211) * qJD(1), 0, 0; 0, 0, 0 (t196 * t200 * t201 + (t198 * t199 * t201 + t195 * t203) * t197) * qJD(3), 0, 0;];
JgD_rot  = t1;
