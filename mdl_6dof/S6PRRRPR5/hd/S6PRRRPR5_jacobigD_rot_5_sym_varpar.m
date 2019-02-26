% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPR5_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:12:56
% EndTime: 2019-02-26 20:12:56
% DurationCPUTime: 0.06s
% Computational Cost: add. (24->17), mult. (88->43), div. (0->0), fcn. (92->10), ass. (0->24)
t193 = sin(pkin(7));
t194 = sin(pkin(6));
t212 = t194 * t193;
t196 = cos(pkin(7));
t200 = cos(qJ(3));
t211 = t196 * t200;
t197 = cos(pkin(6));
t199 = sin(qJ(2));
t210 = t197 * t199;
t201 = cos(qJ(2));
t209 = t197 * t201;
t198 = sin(qJ(3));
t208 = t198 * t201;
t207 = t199 * t200;
t206 = qJD(2) * t198;
t192 = sin(pkin(12));
t195 = cos(pkin(12));
t205 = -t192 * t199 + t195 * t209;
t204 = t192 * t201 + t195 * t210;
t203 = -t192 * t209 - t195 * t199;
t202 = t192 * t210 - t195 * t201;
t191 = t202 * qJD(2);
t190 = t204 * qJD(2);
t1 = [0, 0, -t191 * t193, -t191 * t211 + t203 * t206 + (-t202 * t200 + (t192 * t212 + t203 * t196) * t198) * qJD(3), 0, 0; 0, 0, t190 * t193, t190 * t211 + t205 * t206 + (t204 * t200 + (-t195 * t212 + t205 * t196) * t198) * qJD(3), 0, 0; 0, 0, qJD(2) * t199 * t212, t197 * t193 * qJD(3) * t198 + ((t196 * t208 + t207) * qJD(3) + (t196 * t207 + t208) * qJD(2)) * t194, 0, 0;];
JgD_rot  = t1;
