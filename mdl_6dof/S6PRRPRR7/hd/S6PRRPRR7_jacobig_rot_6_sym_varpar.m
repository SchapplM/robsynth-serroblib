% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRR7_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:37
% EndTime: 2019-02-26 20:07:37
% DurationCPUTime: 0.08s
% Computational Cost: add. (14->9), mult. (39->20), div. (0->0), fcn. (63->8), ass. (0->16)
t104 = sin(pkin(11));
t105 = sin(pkin(6));
t115 = t104 * t105;
t106 = cos(pkin(11));
t114 = t106 * t105;
t107 = cos(pkin(6));
t109 = sin(qJ(2));
t113 = t107 * t109;
t111 = cos(qJ(2));
t112 = t107 * t111;
t110 = cos(qJ(3));
t108 = sin(qJ(3));
t103 = t105 * t109 * t110 + t107 * t108;
t102 = (-t104 * t113 + t106 * t111) * t110 + t108 * t115;
t101 = (t104 * t111 + t106 * t113) * t110 - t108 * t114;
t1 = [0, t115, t104 * t112 + t106 * t109, 0, t102, t102; 0, -t114, t104 * t109 - t106 * t112, 0, t101, t101; 0, t107, -t105 * t111, 0, t103, t103;];
Jg_rot  = t1;
