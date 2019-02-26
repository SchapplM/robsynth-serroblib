% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP12
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRP12_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_jacobigD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:15:32
% EndTime: 2019-02-26 22:15:32
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->16), mult. (78->39), div. (0->0), fcn. (80->8), ass. (0->22)
t188 = sin(pkin(6));
t191 = sin(qJ(2));
t208 = t188 * t191;
t192 = sin(qJ(1));
t207 = t188 * t192;
t195 = cos(qJ(1));
t206 = t188 * t195;
t205 = t191 * t192;
t204 = t191 * t195;
t194 = cos(qJ(2));
t203 = t192 * t194;
t202 = t195 * t194;
t201 = qJD(1) * t188;
t193 = cos(qJ(3));
t200 = qJD(2) * t193;
t189 = cos(pkin(6));
t199 = t189 * t202 - t205;
t198 = t189 * t203 + t204;
t197 = t189 * t204 + t203;
t196 = -t189 * t205 + t202;
t190 = sin(qJ(3));
t1 = [0, t195 * t201, t199 * qJD(1) + t196 * qJD(2), 0 (-t196 * t190 + t193 * t207) * qJD(3) - t198 * t200 + (t190 * t206 - t197 * t193) * qJD(1), 0; 0, t192 * t201, t198 * qJD(1) + t197 * qJD(2), 0 (-t197 * t190 - t193 * t206) * qJD(3) + t199 * t200 + (t190 * t207 + t196 * t193) * qJD(1), 0; 0, 0, qJD(2) * t208, 0, t188 * t194 * t200 + (t189 * t193 - t190 * t208) * qJD(3), 0;];
JgD_rot  = t1;
