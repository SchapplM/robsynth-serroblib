% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRPRR11_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:31
% EndTime: 2019-02-26 20:54:32
% DurationCPUTime: 0.06s
% Computational Cost: add. (24->18), mult. (92->42), div. (0->0), fcn. (96->10), ass. (0->26)
t193 = sin(pkin(7));
t194 = sin(pkin(6));
t214 = t194 * t193;
t196 = cos(pkin(7));
t200 = cos(qJ(3));
t213 = t196 * t200;
t192 = sin(pkin(12));
t199 = sin(qJ(1));
t212 = t199 * t192;
t195 = cos(pkin(12));
t211 = t199 * t195;
t201 = cos(qJ(1));
t210 = t201 * t192;
t209 = t201 * t195;
t208 = t199 * t214;
t207 = t201 * t214;
t206 = qJD(1) * t194 * t196;
t197 = cos(pkin(6));
t205 = t197 * t209 - t212;
t204 = -t197 * t211 - t210;
t203 = t197 * t210 + t211;
t202 = -t197 * t212 + t209;
t198 = sin(qJ(3));
t191 = t204 * qJD(1);
t190 = t205 * qJD(1);
t1 = [0, 0, t190 * t193 + t201 * t206, 0, t190 * t213 + (t202 * t200 + (t204 * t196 + t208) * t198) * qJD(3) + (-t203 * t198 - t200 * t207) * qJD(1), 0; 0, 0, -t191 * t193 + t199 * t206, 0, -t191 * t213 + (t203 * t200 + (t205 * t196 - t207) * t198) * qJD(3) + (t202 * t198 - t200 * t208) * qJD(1), 0; 0, 0, 0, 0 (t193 * t197 * t198 + (t195 * t196 * t198 + t192 * t200) * t194) * qJD(3), 0;];
JgD_rot  = t1;
