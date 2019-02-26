% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRP3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:31:18
% EndTime: 2019-02-26 20:31:18
% DurationCPUTime: 0.03s
% Computational Cost: add. (40->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
t78 = sin(qJ(5));
t79 = sin(qJ(4));
t85 = t79 * t78;
t80 = cos(qJ(5));
t84 = t79 * t80;
t81 = cos(qJ(4));
t83 = t81 * t78;
t82 = t81 * t80;
t77 = qJ(1) + pkin(9);
t76 = cos(t77);
t75 = sin(t77);
t74 = -t75 * t78 + t76 * t84;
t73 = t75 * t80 + t76 * t85;
t72 = t75 * t84 + t76 * t78;
t71 = t75 * t85 - t76 * t80;
t1 = [t74, 0, 0, t75 * t82, -t71, 0; t72, 0, 0, -t76 * t82, t73, 0; 0, 0, 0, -t84, -t83, 0; -t76 * t81, 0, 0, t75 * t79, 0, 0; -t75 * t81, 0, 0, -t76 * t79, 0, 0; 0, 0, 0, t81, 0, 0; t73, 0, 0, t75 * t83, t72, 0; t71, 0, 0, -t76 * t83, -t74, 0; 0, 0, 0, -t85, t82, 0;];
JR_rot  = t1;
