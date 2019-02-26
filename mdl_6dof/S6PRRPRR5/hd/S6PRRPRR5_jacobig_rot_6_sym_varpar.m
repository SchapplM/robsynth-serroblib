% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRR5_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:06:27
% EndTime: 2019-02-26 20:06:28
% DurationCPUTime: 0.08s
% Computational Cost: add. (14->9), mult. (39->20), div. (0->0), fcn. (63->8), ass. (0->16)
t106 = sin(pkin(11));
t107 = sin(pkin(6));
t117 = t106 * t107;
t108 = cos(pkin(11));
t116 = t108 * t107;
t109 = cos(pkin(6));
t111 = sin(qJ(2));
t115 = t109 * t111;
t113 = cos(qJ(2));
t114 = t109 * t113;
t112 = cos(qJ(3));
t110 = sin(qJ(3));
t105 = t107 * t111 * t110 - t109 * t112;
t104 = (-t106 * t115 + t108 * t113) * t110 - t112 * t117;
t103 = (t106 * t113 + t108 * t115) * t110 + t112 * t116;
t1 = [0, t117, t106 * t114 + t108 * t111, 0, t104, t104; 0, -t116, t106 * t111 - t108 * t114, 0, t103, t103; 0, t109, -t107 * t113, 0, t105, t105;];
Jg_rot  = t1;
