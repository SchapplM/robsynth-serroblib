% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:06
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPPRR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:06:01
% EndTime: 2019-02-22 10:06:01
% DurationCPUTime: 0.04s
% Computational Cost: add. (38->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
t67 = sin(qJ(6));
t68 = sin(qJ(5));
t74 = t68 * t67;
t69 = cos(qJ(6));
t73 = t68 * t69;
t70 = cos(qJ(5));
t72 = t70 * t67;
t71 = t70 * t69;
t66 = qJ(1) + pkin(9);
t65 = cos(t66);
t64 = sin(t66);
t63 = -t64 * t67 + t65 * t73;
t62 = -t64 * t69 - t65 * t74;
t61 = -t64 * t73 - t65 * t67;
t60 = t64 * t74 - t65 * t69;
t1 = [t61, 0, 0, 0, t65 * t71, t62; t63, 0, 0, 0, t64 * t71, -t60; 0, 0, 0, 0, -t73, -t72; t60, 0, 0, 0, -t65 * t72, -t63; t62, 0, 0, 0, -t64 * t72, t61; 0, 0, 0, 0, t74, -t71; t64 * t70, 0, 0, 0, t65 * t68, 0; -t65 * t70, 0, 0, 0, t64 * t68, 0; 0, 0, 0, 0, t70, 0;];
JR_rot  = t1;
