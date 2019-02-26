% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR9_jacobig_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobig_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobig_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:27
% EndTime: 2019-02-26 22:20:27
% DurationCPUTime: 0.02s
% Computational Cost: add. (6->6), mult. (17->14), div. (0->0), fcn. (28->8), ass. (0->12)
t105 = cos(pkin(6));
t108 = cos(qJ(2));
t112 = t105 * t108;
t103 = sin(pkin(6));
t107 = sin(qJ(1));
t111 = t107 * t103;
t109 = cos(qJ(1));
t110 = t109 * t103;
t106 = sin(qJ(2));
t104 = cos(pkin(7));
t102 = sin(pkin(7));
t1 = [0, t111 -(-t109 * t106 - t107 * t112) * t102 + t104 * t111, 0, 0, 0; 0, -t110 -(-t107 * t106 + t109 * t112) * t102 - t104 * t110, 0, 0, 0; 1, t105, -t103 * t108 * t102 + t105 * t104, 0, 0, 0;];
Jg_rot  = t1;
